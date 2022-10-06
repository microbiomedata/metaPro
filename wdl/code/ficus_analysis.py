import pandas as pd
import warnings
import sys

warnings.filterwarnings("ignore")

import gffpandas.gffpandas as gffpd
import os

pd.set_option("precision", 20)

__author__ = "Anubhav <anubhav@pnnl.gov>"


class DataOutputtable:
    """
    Created based of data manipulation performed using MSSQLsever queries by
    Purvine, Samuel O <Samuel.Purvine@pnnl.gov>
    """

    def __init__(
        self,
        gff_file,
        resultant_file,
        fasta_txt_file,
        qvalue_threshold,
        dataset_id,
        genome_name,
        dataset_name
    ):

        self.dataset_id = dataset_id
        self.genome_name = genome_name
        self.QValue_threshold = float(qvalue_threshold)
        # read annotation file.
        self.gff_file = gff_file
        self.annotation = gffpd.read_gff3(gff_file)
        self.dataset_name = dataset_name
        self.resultant_df = pd.read_csv(
            resultant_file, sep="\t", float_precision="round_trip"
        )

        col_of_interest = ["ProteinName", "Sequence"]
        self.fasta_df = pd.read_csv(fasta_txt_file, sep="\t")[col_of_interest]
        cols_to_rename = {"ProteinName": "Protein"}
        self.fasta_df.rename(columns=cols_to_rename, inplace=True)
        self.fasta_df["Index"] = self.fasta_df.index + 1

        # saved to served multiple calls
        self.query_8_result = None
        self.peptide_report = None

        # building log excel files:
        file_loc = os.path.join(*self.gff_file.split("/")[:-1])
        # self.writer = pd.ExcelWriter(f"/{file_loc}/{DATA_LOG}.xlsx", engine='xlsxwriter')

        # log dataframe.
        indexes = [f"query_{idx}" for idx in range(23)]
        indexes.extend(["Peptide_Report", "Protein_Report", "qc_metrics"])
        cols = ["#row", "#col", "col_name"]
        self.log_df = pd.DataFrame(columns=cols, index=indexes)

    def _DataOutputtable_xlsx(self, df, sheet_name):

        self.log_df.at[sheet_name, "#row"] = df.shape[0]
        self.log_df.at[sheet_name, "#col"] = df.shape[1]
        self.log_df.at[sheet_name, "col_name"] = ", ".join(df.columns.tolist())
        # df.to_excel(self.writer, sheet_name=sheet_name, index=False)

    def FiltPeps(self, df):
        qvalue_filtered_df = df[df["QValue"] <= self.QValue_threshold]
        non_decoy_filtered_df = qvalue_filtered_df[
            ~qvalue_filtered_df["Protein"].str.startswith("XXX")
        ]
        non_decoy_filtered_df["PeptideSequence"] = non_decoy_filtered_df[
            "Peptide"
        ].str.extract(r"\.(.*)\.")
        return non_decoy_filtered_df

    def get_FiltPeps_gen(self, df):
        FiltPeps_df = self.FiltPeps(df)
        return FiltPeps_df.drop_duplicates().sort_values(by=["PeptideSequence"])

    def get_FiltPeps(self, df):
        FiltPeps_df = self.FiltPeps(df)
        return (
            FiltPeps_df[["PeptideSequence", "Protein"]]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"])
        )

    def query_2(self):
        """
        counting distinct peptides per protein using 5% FDR filter and removing Decoy proteins
        :param query_1_df:
        :return: dataframe with col: 'Protein'  'PeptideSequence_Count'
        """
        temp_df = self.resultant_df.copy()
        FiltPeps_df = self.get_FiltPeps(temp_df)
        grouped_by_protein = FiltPeps_df.groupby("Protein")
        counted_peptides_df = (
            grouped_by_protein["PeptideSequence"]
            .count()
            .reset_index(name="PeptideSequence_Count")
        )

        self._DataOutputtable_xlsx(counted_peptides_df, "query_2")
        return counted_peptides_df

    def query_3(self):
        """

        :param query_1_df:
        :return: dataframe with col: 'PeptideSequence'  'Protein_Count'
        """
        temp_df = self.resultant_df.copy()
        FiltPeps_df = self.get_FiltPeps(temp_df)
        grouped_by_peptide = FiltPeps_df.groupby("PeptideSequence")
        counted_protein_df = (
            grouped_by_peptide["Protein"].count().reset_index(name="Protein_Count")
        )

        self._DataOutputtable_xlsx(counted_protein_df, "query_3")
        return counted_protein_df

    def query_1(self):
        """
        selecting distinct peptide and protein pairs from 5% FDR non-decoy results.
        :param resultant_df:
        :return:
        """
        temp_df = self.resultant_df.copy()
        FiltPeps_df = self.get_FiltPeps(temp_df)

        self._DataOutputtable_xlsx(
            FiltPeps_df.drop_duplicates().sort_values(by=["PeptideSequence"]), "query_1"
        )
        return FiltPeps_df.drop_duplicates().sort_values(by=["PeptideSequence"])

    def query_4(self, table_1, table_2, table_3):

        # inner join on 'Protein'
        merge_1 = table_1.merge(table_2, how="inner", on="Protein")
        # inner join on 'Protein'
        merge_2 = merge_1.merge(self.fasta_df, how="inner", on="Protein")
        # inner join on 'PeptideSequence'
        merge_3 = merge_2.merge(table_3, how="inner", on="PeptideSequence")

        col_of_interest = [
            "Protein",
            "PeptideSequence",
            "PeptideSequence_Count",
            "Index",
            "Protein_Count",
        ]

        self._DataOutputtable_xlsx(
            merge_3[col_of_interest].sort_values(by=["PeptideSequence"]), "query_4"
        )
        return merge_3[col_of_interest].sort_values(by=["PeptideSequence"])

    def get_MaxPepCounts(self, table_4):
        peps_morethan_1_protein = table_4[table_4["Protein_Count"] > 1]
        peps_morethan_1_protein[
            "MaxPeptideSequenceCount"
        ] = peps_morethan_1_protein.groupby(["PeptideSequence"])[
            "PeptideSequence_Count"
        ].transform(
            max
        )
        return peps_morethan_1_protein[
            ["PeptideSequence", "MaxPeptideSequenceCount"]
        ].drop_duplicates()

    def get_MaxPepCounts_table_4_joined(self, table_4):
        MaxPepCounts_df = self.get_MaxPepCounts(table_4)
        joined = MaxPepCounts_df.merge(
            table_4,
            how="inner",
            left_on=["PeptideSequence", "MaxPeptideSequenceCount"],
            right_on=["PeptideSequence", "PeptideSequence_Count"],
        )
        joined["CountOfPeptideCounts"] = joined.groupby(["PeptideSequence"])[
            "PeptideSequence_Count"
        ].transform("count")
        return joined

    def get_BestSingleProteins(self, table_4):

        MaxPepCounts_table_4_joined_df = self.get_MaxPepCounts_table_4_joined(table_4)
        filtered = MaxPepCounts_table_4_joined_df[
            MaxPepCounts_table_4_joined_df["CountOfPeptideCounts"] == 1
        ]
        filtered["MaxPeptideSequenceCount"] = filtered.groupby(["PeptideSequence"])[
            "PeptideSequence_Count"
        ].transform(max)
        return filtered[
            ["PeptideSequence", "CountOfPeptideCounts", "MaxPeptideSequenceCount"]
        ]

    def get_BestIndexedProtein(self, table_4):
        MaxPepCounts_table_4_joined_df = self.get_MaxPepCounts_table_4_joined(table_4)
        filtered = MaxPepCounts_table_4_joined_df[
            MaxPepCounts_table_4_joined_df["CountOfPeptideCounts"] > 1
        ]
        filtered["MinIndexValue"] = filtered.groupby(["PeptideSequence"])[
            "Index"
        ].transform(min)
        return filtered[["PeptideSequence", "CountOfPeptideCounts", "MinIndexValue"]]

    def query_5(self, table_4):
        peps_with_1_protein = table_4[table_4["Protein_Count"] == 1]
        peps_with_1_protein.rename(columns={"Protein": "BestProtein"}, inplace=True)

        self._DataOutputtable_xlsx(
            peps_with_1_protein[["PeptideSequence", "BestProtein"]], "query_5"
        )
        return peps_with_1_protein[["PeptideSequence", "BestProtein"]]

    def query_6(self, table_4):
        BestSingleProteins_df = self.get_BestSingleProteins(table_4)
        joined = BestSingleProteins_df.merge(
            table_4,
            how="inner",
            left_on=["PeptideSequence", "MaxPeptideSequenceCount"],
            right_on=["PeptideSequence", "PeptideSequence_Count"],
        )
        joined.rename(columns={"Protein": "BestProtein"}, inplace=True)

        self._DataOutputtable_xlsx(
            joined[["PeptideSequence", "BestProtein"]], "query_6"
        )
        return joined[["PeptideSequence", "BestProtein"]]

    def query_7(self, table_4):
        BestIndexedProtein_df = self.get_BestIndexedProtein(table_4)
        merged = BestIndexedProtein_df.merge(
            table_4,
            how="inner",
            left_on=["PeptideSequence", "MinIndexValue"],
            right_on=["PeptideSequence", "Index"],
        )
        merged.rename(columns={"Protein": "BestProtein"}, inplace=True)

        self._DataOutputtable_xlsx(
            merged[["PeptideSequence", "BestProtein"]].drop_duplicates(), "query_7"
        )
        return merged[["PeptideSequence", "BestProtein"]].drop_duplicates()

    def get_best_protein_associated_with_peptide(self, table_4):
        # Take Union of 3 dataframes!
        table_5 = self.query_5(table_4)
        table_6 = self.query_6(table_4)
        table_7 = self.query_7(table_4)
        return pd.concat([table_5, table_6, table_7])

    def parse_MSGFjobs_MASIC_resultant(self):
        table_1 = self.query_1()
        table_2 = self.query_2()
        table_3 = self.query_3()
        table_4 = self.query_4(table_1, table_2, table_3)
        table_5_6_7_unioned = self.get_best_protein_associated_with_peptide(table_4)
        return table_5_6_7_unioned

    def query_0(self):
        """

        :param gff_file: JGI annotation file.
        :return: dataframe with col: ['protein', 'Product', 'pfam', 'ko', 'ec_number', 'cog']
        """
        filtered_df = self.annotation.filter_feature_of_type(["CDS"])
        col_of_interest = ["ID", "product", "pfam", "ko", "ec_number", "cog"]
        attributeTag_df = filtered_df.attributes_to_columns()[col_of_interest]
        cols_to_rename = {"ID": "Protein", "product": "Product"}
        attributeTag_df.rename(columns=cols_to_rename, inplace=True)
        # cleaning up columns
        # TODO: remove KO: EC: strings
        attributeTag_df["ko"] = attributeTag_df["ko"]  # .str.replace('KO:', '')
        attributeTag_df["ec_number"] = attributeTag_df[
            "ec_number"
        ]  # .str.replace('EC:', '')

        self._DataOutputtable_xlsx(attributeTag_df, "query_0")
        return attributeTag_df

    def query_8(self):
        # TODO:remove all tables except table_5_6_7_unioned
        # 1. get best_protein_associated_with_peptide
        #         table_1,table_2,table_3,table_4,
        table_5_6_7_unioned = self.parse_MSGFjobs_MASIC_resultant()
        # 2. get annotations
        AnnotationSplitout_df = self.query_0()

        merged = table_5_6_7_unioned.merge(
            AnnotationSplitout_df,
            how="left",
            left_on=["BestProtein"],
            right_on=["Protein"],
        )
        col_of_interest = [
            "PeptideSequence",
            "BestProtein",
            "Product",
            "ec_number",
            "pfam",
            "ko",
            "cog",
        ]

        cols_to_rename = {"ec_number": "EC_Number", "ko": "KO", "cog": "COG"}
        PeptideBestProteinAnnotated = merged[col_of_interest]
        PeptideBestProteinAnnotated.rename(columns=cols_to_rename, inplace=True)

        self._DataOutputtable_xlsx(PeptideBestProteinAnnotated, "query_8")
        return PeptideBestProteinAnnotated

    def get_FwdPeptideSequences(self):
        FiltPeps_df = self.get_FiltPeps(self.resultant_df)
        FiltPeps_df["GeneCount"] = FiltPeps_df.groupby(["PeptideSequence"])[
            "Protein"
        ].transform("nunique")
        return FiltPeps_df[["PeptideSequence", "GeneCount"]]

    def query_9(self):
        FwdPeptideSequences_df = self.get_FwdPeptideSequences()

        temp_df = self.resultant_df.copy()
        # temp_df['Peptide'] = temp_df['Peptide'].str.extract(r'\.(.*)\.')

        FiltPassPeptideAbundanceData = self.get_FiltPeps_gen(temp_df)[
            ["PeptideSequence", "Protein"]
        ]

        merged = FwdPeptideSequences_df.merge(
            FiltPassPeptideAbundanceData,
            how="inner",
            left_on=["PeptideSequence"],
            right_on=["PeptideSequence"],
        )

        self._DataOutputtable_xlsx(
            merged[["PeptideSequence", "GeneCount", "Protein"]]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"]),
            "query_9",
        )
        return (
            merged[["PeptideSequence", "GeneCount", "Protein"]]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"])
        )

    def build_annotation_str(self, Protein, Product, pfam, ko, ec_number, cog):

        if not Protein.startswith("Contaminant"):

            annotation_str = f"gene_name={Protein}"
            if Product and (Product != None) and (Product != "nan"):
                annotation_str = annotation_str + f";product={Product}"
            if pfam and (pfam != None) and (pfam != "nan"):
                annotation_str = annotation_str + f";pfam={pfam}"
            if ko and (ko != None) and (ko != "nan"):
                annotation_str = annotation_str + f";ko={ko}"
            if ec_number and (ec_number != None) and (ec_number != "nan"):
                annotation_str = annotation_str + f";ec_number={ec_number}"
            if cog and (cog != None) and (cog != "nan"):
                annotation_str = annotation_str + f";cog={cog}"
        else:
            # the contaminants aren't present in the JGI provided annotation files!
            # contaminants are added from PNNL.
            annotation_str = Protein
        return annotation_str

    def query_10(self):
        table_9 = self.query_9()
        AnnotationSplitout_df = self.query_0()
        merged = table_9.merge(
            AnnotationSplitout_df, how="left", left_on=["Protein"], right_on=["Protein"]
        )
        merged["annotation"] = merged.apply(
            lambda x: self.build_annotation_str(
                x.Protein, x.Product, x.pfam, x.ko, x.ec_number, x.cog
            ),
            axis=1,
        )

        self._DataOutputtable_xlsx(
            merged[
                ["PeptideSequence", "GeneCount", "Protein", "annotation"]
            ].sort_values(by=["PeptideSequence"]),
            "query_10",
        )
        return merged[
            ["PeptideSequence", "GeneCount", "Protein", "annotation"]
        ].sort_values(by=["PeptideSequence"])

    def query_11(self):
        table_10 = self.query_10()

        # prerare protein list
        table_10["FullGeneList"] = (
            table_10.sort_values("Protein")
            .groupby(["PeptideSequence", "GeneCount"])["Protein"]
            .transform(lambda x: ", ".join(x))
        )

        # sort them!
        get_list_func = lambda x: x.split(", ") if len(x.split(", ")) > 1 else x
        f = (
            lambda x: ", ".join(sorted(get_list_func(x)))
            if isinstance(get_list_func(x), list)
            else get_list_func(x)
        )
        table_10["FullGeneList"] = table_10["FullGeneList"].apply(f)

        table_10["AnnotationList"] = (
            table_10.sort_values("annotation")
            .groupby(["PeptideSequence", "GeneCount"])["annotation"]
            .transform(lambda x: " | ".join(x))
        )

        self._DataOutputtable_xlsx(
            table_10[["PeptideSequence", "GeneCount", "FullGeneList", "AnnotationList"]]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"]),
            "query_11",
        )
        return (
            table_10[["PeptideSequence", "GeneCount", "FullGeneList", "AnnotationList"]]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"])
        )

    def query_12(self):

        temp_df = self.resultant_df.copy()
        FiltPassPeptideAbundanceData = self.get_FiltPeps_gen(temp_df)[
            ["SpecIndex", "PeptideSequence", "StatMomentsArea"]
        ].drop_duplicates(
            subset=["SpecIndex", "PeptideSequence", "StatMomentsArea"]
        )  # .sort_values(by=[''])
        FiltPassPeptideAbundanceData[
            "SpectralCount"
        ] = FiltPassPeptideAbundanceData.groupby(["PeptideSequence"])[
            "SpecIndex"
        ].transform(
            "count"
        )
        FiltPassPeptideAbundanceData[
            "sum(StatMomentsArea)"
        ] = FiltPassPeptideAbundanceData.groupby(["PeptideSequence"])[
            "StatMomentsArea"
        ].transform(
            "sum"
        )

        self._DataOutputtable_xlsx(
            FiltPassPeptideAbundanceData[
                ["PeptideSequence", "SpectralCount", "sum(StatMomentsArea)"]
            ]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"]),
            "query_12",
        )
        return (
            FiltPassPeptideAbundanceData[
                ["PeptideSequence", "SpectralCount", "sum(StatMomentsArea)"]
            ]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"])
        )

    def query_13(self):
        temp_df = self.resultant_df.copy()
        # temp_df['Dataset_x'] = 'default value'
        filtered = self.get_FiltPeps_gen(temp_df)[
            ["PeptideSequence", "QValue"]
            # ["Dataset_x", "PeptideSequence", "QValue"]
        ]
        groupby_columns = ["PeptideSequence"]
        # groupby_columns = ["Dataset_x", "PeptideSequence"]
        filtered["min(QValue)"] = filtered.groupby(groupby_columns)["QValue"].transform(
            "min"
        )
        PeptideGeneListsANDAnnotationList_df = self.query_11()

        # inner join on 'PeptideSequence'
        merged_1 = filtered.merge(
            PeptideGeneListsANDAnnotationList_df,
            how="inner",
            left_on=["PeptideSequence"],
            right_on=["PeptideSequence"],
        )

        PeptideBestProteinAnnotated_df = self.query_8()
        # saved for multiple call.
        self.query_8_result = PeptideBestProteinAnnotated_df.copy()

        # left outer join 'PeptideSequence'
        merged_2 = merged_1.merge(
            PeptideBestProteinAnnotated_df, how="left", on=["PeptideSequence"]
        )

        PeptideAbundanceInfo_df = self.query_12()
        # inner join on 'PeptideSequence'
        merged_3 = merged_2.merge(
            PeptideAbundanceInfo_df,
            how="inner",
            left_on=["PeptideSequence"],
            right_on=["PeptideSequence"],
        )

        cols_to_rename = {
            # "Dataset_x": "DatasetName",
            "sum(StatMomentsArea)": "sum(MASICAbundance)",
        }
        merged_3.rename(columns=cols_to_rename, inplace=True)

        merged_3["DatasetName"] = self.dataset_name

        self._DataOutputtable_xlsx(
            merged_3[
                [
                    "DatasetName",
                    "PeptideSequence",
                    "BestProtein",
                    "Product",
                    "EC_Number",
                    "pfam",
                    "KO",
                    "COG",
                    "GeneCount",
                    "FullGeneList",
                    "AnnotationList",
                    "min(QValue)",
                    "SpectralCount",
                    "sum(MASICAbundance)",
                ]
            ]
            .drop_duplicates()
            .sort_values(by=["DatasetName"]),
            "query_13",
        )
        return (
            merged_3[
                [
                    "DatasetName",
                    "PeptideSequence",
                    "BestProtein",
                    "Product",
                    "EC_Number",
                    "pfam",
                    "KO",
                    "COG",
                    "GeneCount",
                    "FullGeneList",
                    "AnnotationList",
                    "min(QValue)",
                    "SpectralCount",
                    "sum(MASICAbundance)",
                ]
            ]
            .drop_duplicates()
            .sort_values(by=["DatasetName"])
        )

    def query_14(self):
        Peptide_Report = self.query_13()

        Peptide_Report["UniquePeptideCount"] = Peptide_Report.groupby(["BestProtein"])[
            "PeptideSequence"
        ].transform("count")
        Peptide_Report["SummedSpectraCounts"] = Peptide_Report.groupby(["BestProtein"])[
            "SpectralCount"
        ].transform("sum")
        Peptide_Report["SummedPeptideMASICAbundances"] = Peptide_Report.groupby(
            ["BestProtein"]
        )["sum(MASICAbundance)"].transform("sum")

        self._DataOutputtable_xlsx(
            Peptide_Report[
                [
                    "BestProtein",
                    "UniquePeptideCount",
                    "SummedSpectraCounts",
                    "SummedPeptideMASICAbundances",
                ]
            ]
            .drop_duplicates()
            .sort_values(by=["BestProtein"]),
            "query_14",
        )
        return (
            Peptide_Report[
                [
                    "BestProtein",
                    "UniquePeptideCount",
                    "SummedSpectraCounts",
                    "SummedPeptideMASICAbundances",
                ]
            ]
            .drop_duplicates()
            .sort_values(by=["BestProtein"])
        )

    def query_15(self):

        PeptideBestProteinAnnotated_df = self.query_8_result  # self.query_8()
        temp_df = self.resultant_df.copy()
        temp_df["PeptideSequence"] = temp_df["Peptide"].str.extract(r"\.(.*)\.")
        merged = PeptideBestProteinAnnotated_df.merge(
            temp_df,
            how="inner",
            left_on=["PeptideSequence"],
            right_on=["PeptideSequence"],
        )
        non_decoy_filtered_df = merged[~merged["Protein"].str.startswith("XXX")]
        self._DataOutputtable_xlsx(
            non_decoy_filtered_df[["BestProtein", "Protein"]]
            .drop_duplicates(subset=["BestProtein", "Protein"])
            .sort_values(by=["BestProtein"]),
            "query_15",
        )
        return (
            non_decoy_filtered_df[["BestProtein", "Protein"]]
            .drop_duplicates(subset=["BestProtein", "Protein"])
            .sort_values(by=["BestProtein"])
        )

    def query_16(self):
        BestProteinGeneLists = self.query_15()  # .drop_duplicates(subset=['Protein'])

        # prepare protein list
        BestProteinGeneLists["FullGeneList"] = BestProteinGeneLists.groupby(
            ["BestProtein"]
        )["Protein"].transform(lambda x: ", ".join(x))
        # sort them!
        get_list_func = lambda x: x.split(", ") if len(x.split(", ")) > 1 else x
        f = (
            lambda x: ", ".join(sorted(get_list_func(x)))
            if isinstance(get_list_func(x), list)
            else get_list_func(x)
        )
        BestProteinGeneLists["FullGeneList"] = BestProteinGeneLists[
            "FullGeneList"
        ].apply(f)

        BestProteinGeneLists["GeneCount"] = BestProteinGeneLists.groupby(
            ["BestProtein"]
        )["Protein"].transform("count")
        BestProteinGeneLists.sort_values(["FullGeneList"], ascending=True, inplace=True)

        self._DataOutputtable_xlsx(
            BestProteinGeneLists[["BestProtein", "GeneCount", "FullGeneList"]]
            .drop_duplicates()
            .sort_values(by=["BestProtein"]),
            "query_16",
        )
        return (
            BestProteinGeneLists[["BestProtein", "GeneCount", "FullGeneList"]]
            .drop_duplicates()
            .sort_values(by=["BestProtein"])
        )

    def mod_build_annotation_str(
        self, BestProtein, Protein, Product, pfam, ko, ec_number, cog
    ):

        if Product != Product:  # pd.isna(Product):
            # handle <class 'float'> nan
            return BestProtein
        else:

            annotation_str = f"gene_name={Protein}" + f";product={Product}"
            if pfam and (pfam != None) and (pfam != "nan"):
                annotation_str = annotation_str + f";pfam={pfam}"
            if ko and (ko != None) and (ko != "nan"):
                annotation_str = annotation_str + f";ko={ko}"
            if ec_number and (ec_number != None) and (ec_number != "nan"):
                annotation_str = annotation_str + f";ec_number={ec_number}"
            if cog and (cog != None) and (cog != "nan"):
                annotation_str = annotation_str + f";cog={cog}"
            return annotation_str

    def query_17(self):

        AnnotationSplitout_df = self.query_0()

        table_15 = self.query_15()  # have bestProtein
        ProteinBestProteinsAnnotated = table_15.merge(
            AnnotationSplitout_df, how="left", left_on=["Protein"], right_on=["Protein"]
        )
        #         display(ProteinBestProteinsAnnotated)
        ProteinBestProteinsAnnotated["annotation"] = ProteinBestProteinsAnnotated.apply(
            lambda x: self.mod_build_annotation_str(
                x.BestProtein, x.Protein, x.Product, x.pfam, x.ko, x.ec_number, x.cog
            ),
            axis=1,
        )

        # display(ProteinBestProteinsAnnotated)
        self._DataOutputtable_xlsx(
            ProteinBestProteinsAnnotated[["BestProtein", "Protein", "annotation"]]
            .drop_duplicates()
            .sort_values(by=["BestProtein"]),
            "query_17",
        )

        return (
            ProteinBestProteinsAnnotated[["BestProtein", "Protein", "annotation"]]
            .drop_duplicates()
            .sort_values(by=["BestProtein"])
        )

    def query_18(self):

        BestProteinAnnotationLists = self.query_17()

        BestProteinAnnotationLists[
            "AnnotationList"
        ] = BestProteinAnnotationLists.groupby(["BestProtein"])["annotation"].transform(
            lambda x: " | ".join(x)
        )
        self._DataOutputtable_xlsx(
            BestProteinAnnotationLists[["BestProtein", "AnnotationList"]]
            .drop_duplicates()
            .sort_values(by=["BestProtein"]),
            "query_18",
        )
        return (
            BestProteinAnnotationLists[["BestProtein", "AnnotationList"]]
            .drop_duplicates()
            .sort_values(by=["BestProtein"])
        )

    def query_19(self):
        Peptide_Report_df = self.peptide_report.copy()[
            ["DatasetName", "BestProtein", "Product", "EC_Number", "pfam", "KO", "COG"]
        ]
        BestProteinGeneLists = self.query_16()
        merged_1 = Peptide_Report_df.merge(
            BestProteinGeneLists,
            how="inner",
            left_on=["BestProtein"],
            right_on=["BestProtein"],
        )

        BestProteinAnnotationLists = self.query_18()
        merged_2 = merged_1.merge(
            BestProteinAnnotationLists,
            how="inner",
            left_on=["BestProtein"],
            right_on=["BestProtein"],
        )

        PeptideAbundanceInfo_df = self.query_14()
        merged_3 = merged_2.merge(
            PeptideAbundanceInfo_df,
            how="inner",
            left_on=["BestProtein"],
            right_on=["BestProtein"],
        )

        distinct_colm = [
            "DatasetName",
            "BestProtein",
            "Product",
            "EC_Number",
            "pfam",
            "KO",
            "COG",
            "GeneCount",
            "FullGeneList",
            "AnnotationList",
            "UniquePeptideCount",
            "SummedSpectraCounts",
            "SummedPeptideMASICAbundances",
        ]
        self._DataOutputtable_xlsx(
            merged_3[distinct_colm]
            .drop_duplicates(subset=distinct_colm)
            .sort_values(by=["DatasetName"]),
            "query_19",
        )
        return (
            merged_3[distinct_colm]
            .drop_duplicates(subset=distinct_colm)
            .sort_values(by=["DatasetName"])
        )

    def query_20(self):

        total_protein_count = pd.DataFrame(
            {"total_protein_count": self.fasta_df.Protein.nunique()}, index=[0]
        )
        self._DataOutputtable_xlsx(total_protein_count, "query_20")
        return total_protein_count

    def query_21(self):
        temp_df = self.resultant_df.copy()
        filtered = self.get_FiltPeps_gen(temp_df)[["Scan", "Charge"]]
        PSM_Count = pd.DataFrame(
            {"PSM_Count": len(filtered.groupby(["Scan", "Charge"])["Scan"])}, index=[0]
        )
        self._DataOutputtable_xlsx(PSM_Count, "query_21")
        return PSM_Count

    def query_22(self):
        temp_df = self.resultant_df.copy()
        Total_PSM_Count = pd.DataFrame(
            {
                "Total_PSM_Count": len(
                    temp_df[["Scan", "Charge"]].groupby(["Scan", "Charge"])["Scan"]
                )
            },
            index=[0],
        )
        self._DataOutputtable_xlsx(Total_PSM_Count, "query_22")
        return Total_PSM_Count

    def query_23(self):
        Table12 = pd.DataFrame(
            {
                "PSM_identification_rate": (
                    self.query_21().PSM_Count / self.query_22().Total_PSM_Count
                )
                * 100
            },
            index=[0],
        )
        self._DataOutputtable_xlsx(Table12, "query_23")
        return Table12

    def query_24(self):

        Peptide_Report_df = self.peptide_report.copy()
        unique_peptide_seq_count = (
            Peptide_Report_df.PeptideSequence.drop_duplicates().count()
        )
        BestProtein_count = Peptide_Report_df.BestProtein.drop_duplicates().count()
        mean_peptide_count = (unique_peptide_seq_count / BestProtein_count) * 100

        qc_metrics = pd.DataFrame(
            {
                "total_protein_count": self.query_20().total_protein_count,
                "PSM_count": self.query_21().PSM_Count,
                "PSM_identification_rate": self.query_23().PSM_identification_rate,
                "unique_peptide_seq_count": unique_peptide_seq_count,
                "BestProtein_count": BestProtein_count,
                "mean_peptide_count": mean_peptide_count,
            },
            index=[0],
        )

        self._DataOutputtable_xlsx(qc_metrics, "query_24")
        return qc_metrics

    def create_peptide_report(self):
        Peptide_Report_df = self.query_13()
        Peptide_Report_df.sort_values(["GeneCount"], ascending=False, inplace=True)
        self._DataOutputtable_xlsx(Peptide_Report_df, "Peptide_Report")
        return Peptide_Report_df

    def create_protein_report(self):
        Protein_Report_df = self.query_19()
        Protein_Report_df.sort_values(["GeneCount"], ascending=False, inplace=True)

        self._DataOutputtable_xlsx(Protein_Report_df, "Protein_Report")
        return Protein_Report_df

    def create_qc_metrics(self):
        qc_metrics_df = self.query_24()
        self._DataOutputtable_xlsx(qc_metrics_df, "qc_metrics")
        return qc_metrics_df

    def gen_reports(self):
        print(f"ReportGen start {self.dataset_id} : {self.genome_name}")
        self.peptide_report = self.create_peptide_report()
        protein_report = self.create_protein_report()
        qc_metrics_report = self.create_qc_metrics()
        print(f"ReportGen end {self.dataset_id} : {self.genome_name}")

        # rename to adhere nmdc.schema.json
        cols_to_rename = {
            "PeptideSequence": "peptide_sequence",
            "BestProtein": "best_protein",
            "FullGeneList": "all_proteins",
            "min(QValue)": "min_q_value",
            "SpectralCount": "peptide_spectral_count",
            "sum(MASICAbundance)": "peptide_sum_masic_abundance",
        }
        self.peptide_report.rename(columns=cols_to_rename, inplace=True)

        # rename to adhere nmdc.schema.json
        cols_to_rename = {
            "BestProtein": "best_protein",
            "FullGeneList": "all_proteins",
            "sum(MASICAbundance)": "peptide_sum_masic_abundance",
        }
        protein_report.rename(columns=cols_to_rename, inplace=True)

        # TODO: Change string type to "min_q_value", 'peptide_spectral_count', 'peptide_sum_masic_abundance' int!
        return self.peptide_report, protein_report, qc_metrics_report

if __name__ == "__main__":

    fasta_txt_file = sys.argv[1]
    gff_file = sys.argv[2]
    resultant_file = sys.argv[3]
    dataset_id= sys.argv[4]
    genome_directory= sys.argv[5]
    QVALUE_THRESHOLD= sys.argv[6]
    dataset_name= sys.argv[7]

    print(f"{fasta_txt_file}\n{gff_file}\n{resultant_file}\n{dataset_id}\n{genome_directory}\n{QVALUE_THRESHOLD}\n")

    data_obj = DataOutputtable(
        gff_file,
        resultant_file,
        fasta_txt_file,
        QVALUE_THRESHOLD,
        dataset_id,
        genome_directory,
        dataset_name
    )

    (
        peptide_report,
        protein_report,
        qc_metrics_report,
    ) = data_obj.gen_reports()

    # write dfs to file.
    peptide_report.to_csv(f"{dataset_id}_{genome_directory}_Peptide_Report.tsv", sep="\t")
    protein_report.to_csv(f"{dataset_id}_{genome_directory}_Protein_Report.tsv", sep="\t")
    qc_metrics_report.to_csv( f"{dataset_id}_{genome_directory}_QC_metrics.tsv", sep="\t")

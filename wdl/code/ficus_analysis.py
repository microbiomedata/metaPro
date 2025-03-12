import pandas as pd
import warnings
import sys
import json
import functools
from pathlib import Path
from typing import Tuple, Optional
from datetime import datetime

warnings.filterwarnings("ignore")

import gffpandas.gffpandas as gffpd
import os

pd.set_option("display.precision", 20)

def timestamp_as_string():
    return datetime.now().strftime("%Y%m%d%H%M%S")
    
class DataOutputtable:
    """
    Created based of data manipulation performed using MSSQLsever queries by
    Purvine, Samuel O <Samuel.Purvine@pnnl.gov>
    """

    def __init__(
        self,
        gff_file,
        resultant_file,
        fasta_txt_file, # this can go away
        threshold,
        dataset_id,     # this can go away
        faa_id,         # this can go away
        dataset_name,
        did_split_analysis,
        is_metagenome_free_analysis # reevaluate if this is needed or not since new kaiko writes much better .gffs
    ):

        self.dataset_id = dataset_id
        self.faa_id = faa_id
        self.threshold = float(threshold)
        # read annotation file.
        self.gff_file = gff_file
        self.annotation = gffpd.read_gff3(gff_file)
        self.dataset_name = dataset_name
        self.resultant_df = pd.read_csv(resultant_file, sep="\t", float_precision="round_trip")
        self.did_split_analysis = did_split_analysis
        self.is_metagenome_free_analysis = is_metagenome_free_analysis

        # saved to served multiple calls
        self.query_8_result = None
        self.peptide_report = None

        # # building log excel files:
        # gff_parent_path = Path(gff_file).parent
        # self.writer = pd.ExcelWriter(gff_parent_path / f"query_results_{timestamp_as_string()}.xlsx")
        
    def write_df_excel(func):
        '''
        This function wraps all query functions to write the resultant 10k rows of the returned 
        dataframe to an excel file.
        '''
        @functools.wraps(func)
        def capture(self, *args, **kwargs):
            result_df: pd.DataFrame = func(self, *args, **kwargs)   
            # result_df.iloc[:10000].to_excel(self.writer, sheet_name=func.__name__, index=False)
            return result_df
        return capture

    def FiltPeps(self, df: pd.DataFrame) -> pd.DataFrame:

        qvalue_filtered_df = df[df["MSGFDB_SpecEValue"] <= self.threshold] if self.did_split_analysis else df[df["QValue"] <= self.threshold]
        non_decoy_filtered_df = qvalue_filtered_df[
            ~qvalue_filtered_df["Protein"].str.startswith("XXX")
        ]
        non_decoy_filtered_df["PeptideSequence"] = non_decoy_filtered_df[
            "Peptide"
        ].str.extract(r"\.(.*)\.")
        return non_decoy_filtered_df

    @write_df_excel
    def get_FiltPeps_gen(self, df):
        FiltPeps_df = self.FiltPeps(df)
        return FiltPeps_df.drop_duplicates().sort_values(by=["PeptideSequence"])

    @write_df_excel
    def get_FiltPeps(self, df: pd.DataFrame) -> pd.DataFrame:
        FiltPeps_df = self.FiltPeps(df)
        return (
            FiltPeps_df[["PeptideSequence", "Protein"]]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"])
        )

    @write_df_excel
    def query_2(self) -> pd.DataFrame:
        """
        counting distinct peptides per protein using 5% FDR filter and removing Decoy proteins
        :param query_1_df:
        :return: dataframe with col: 'Protein'  'peptides_per_protein'
        """
        temp_df = self.resultant_df.copy()
        FiltPeps_df = self.get_FiltPeps(temp_df)
        grouped_by_protein = FiltPeps_df.groupby("Protein")
        counted_peptides_df = (
            grouped_by_protein["PeptideSequence"]
            .count()
            .reset_index(name="peptides_per_protein")
        )

        return counted_peptides_df

    @write_df_excel
    def query_3(self) -> pd.DataFrame:
        """

        :param query_1_df:
        :return: dataframe with col: 'PeptideSequence'  'proteins_per_peptide'
        """
        temp_df = self.resultant_df.copy()
        FiltPeps_df = self.get_FiltPeps(temp_df)
        grouped_by_peptide = FiltPeps_df.groupby("PeptideSequence")
        counted_protein_df = (
            grouped_by_peptide["Protein"].count().reset_index(name="proteins_per_peptide")
        )

        return counted_protein_df

    @write_df_excel
    def query_1(self) -> pd.DataFrame:
        """
        selecting distinct peptide and protein pairs from 5% FDR non-decoy results.
        :param resultant_df:
        :return:
        """
        temp_df = self.resultant_df.copy()
        FiltPeps_df = self.get_FiltPeps(temp_df)

        return FiltPeps_df.drop_duplicates().sort_values(by=["PeptideSequence"])

    @write_df_excel
    def query_4(self, table_1: pd.DataFrame, table_2: pd.DataFrame, table_3: pd.DataFrame) -> pd.DataFrame:

        # inner join on 'Protein'
        merge_1 = table_1.merge(table_2, how="inner", on="Protein")

        # inner join on 'PeptideSequence'
        merge_3 = merge_1.merge(table_3, how="inner", on="PeptideSequence")

        col_of_interest = [
            "Protein",
            "PeptideSequence",
            "peptides_per_protein",
            "proteins_per_peptide",
        ]

        return merge_3[col_of_interest].sort_values(by=["PeptideSequence"])

    @write_df_excel
    def get_MaxPepCounts(self, table_4: pd.DataFrame) -> pd.DataFrame:
        peps_morethan_1_protein = table_4[table_4["proteins_per_peptide"] > 1]
        peps_morethan_1_protein[
            "max_peptides_per_protein"
        ] = peps_morethan_1_protein.groupby(["PeptideSequence"])[
            "peptides_per_protein"
        ].transform(
            max
        )
        return peps_morethan_1_protein[
            ["PeptideSequence", "max_peptides_per_protein"]
        ].drop_duplicates()

    @write_df_excel
    def get_MaxPepCounts_table_4_joined(self, table_4: pd.DataFrame) -> pd.DataFrame:
        MaxPepCounts_df = self.get_MaxPepCounts(table_4)
        joined = MaxPepCounts_df.merge(
            table_4,
            how="inner",
            left_on=["PeptideSequence", "max_peptides_per_protein"],
            right_on=["PeptideSequence", "peptides_per_protein"],
        )
        joined["CountOfPeptideCounts"] = joined.groupby(["PeptideSequence"])[
            "peptides_per_protein"
        ].transform("count")
        return joined

    @write_df_excel
    def get_BestIndexedProtein(self, table_4: pd.DataFrame) -> pd.DataFrame:
        MaxPepCounts_table_4_joined_df = self.get_MaxPepCounts_table_4_joined(table_4)
        filtered = MaxPepCounts_table_4_joined_df[
            MaxPepCounts_table_4_joined_df["CountOfPeptideCounts"] > 1
        ]

        return filtered[["PeptideSequence", "CountOfPeptideCounts"]]

    @write_df_excel
    def get_razor_protein_associated_with_peptide(self, table_4: pd.DataFrame) -> pd.DataFrame:

        peps_with_1_protein = table_4[table_4["proteins_per_peptide"] == 1]

        non_degenerate_razor_protein_df = peps_with_1_protein[["PeptideSequence", "Protein"]] 

        peps_with_many_proteins = table_4[table_4["proteins_per_peptide"] > 1]

        non_degenerate_protein_set = set(non_degenerate_razor_protein_df['Protein'].unique())
        NonDegenerate_razor_protein_list_member_int_df = peps_with_many_proteins[peps_with_many_proteins['Protein'].isin(non_degenerate_protein_set)]

        Degenerate_peptides_with_only_unique_protein =  NonDegenerate_razor_protein_list_member_int_df.drop_duplicates(subset=['PeptideSequence'], keep=False, inplace=False)[['PeptideSequence', 'Protein']]

        Degenerate_peptides_with_morethanone_unique_protein_to_remove_df = NonDegenerate_razor_protein_list_member_int_df[NonDegenerate_razor_protein_list_member_int_df['PeptideSequence'].duplicated(keep=False)]

        peptides_found_or_unable_to_match_set = set(Degenerate_peptides_with_morethanone_unique_protein_to_remove_df['PeptideSequence']) \
            | set(Degenerate_peptides_with_only_unique_protein['PeptideSequence']) \
            | set(non_degenerate_razor_protein_df['PeptideSequence'])
        
        max_peptide_per_protein_count = peps_with_many_proteins.copy()
        max_peptide_per_protein_count = max_peptide_per_protein_count[~max_peptide_per_protein_count['PeptideSequence'].isin(peptides_found_or_unable_to_match_set)]
        max_peptide_per_protein_count["max_peptides_per_protein"] = max_peptide_per_protein_count.groupby(["PeptideSequence"])["peptides_per_protein"].transform(max)
        max_peptide_per_protein_count = max_peptide_per_protein_count[max_peptide_per_protein_count["peptides_per_protein"] == max_peptide_per_protein_count["max_peptides_per_protein"]]
        Degenerate_peptides_with_max_nonunique_peptide_count_df = max_peptide_per_protein_count[["PeptideSequence", "Protein"]]

        result = pd.concat([non_degenerate_razor_protein_df, Degenerate_peptides_with_only_unique_protein, Degenerate_peptides_with_max_nonunique_peptide_count_df])
        result = result.rename(columns={"Protein": "RazorProtein"})
        return result

    @write_df_excel
    def parse_MSGFjobs_MASIC_resultant(self) -> pd.DataFrame:
        table_1 = self.query_1()
        table_2 = self.query_2()
        table_3 = self.query_3()
        table_4 = self.query_4(table_1, table_2, table_3)
        peptide_to_razor_protein_df = self.get_razor_protein_associated_with_peptide(table_4)
        return peptide_to_razor_protein_df

    @write_df_excel
    def query_0(self) -> pd.DataFrame:
        """

        :param gff_file: JGI annotation file.
        :return: dataframe with col: ['protein', 'Product', 'pfam', 'ko', 'ec_number', 'cog']
        """

        if self.is_metagenome_free_analysis:
            df = pd.DataFrame()
            col_of_interest = ["ID", "product", "kegg_annotations"]
            df["parsed_attributes"] = self.annotation.df["attributes"].apply(self.parse_attributes)
            attributeTag_df = pd.json_normalize(df["parsed_attributes"]).fillna("None")[col_of_interest]
            attributeTag_df["pfam"] = None
            attributeTag_df["ec_number"] = None
            attributeTag_df["cog"] = None
            cols_to_rename = {"ID": "Protein", "product": "Product", "kegg_annotations": "ko"}
            attributeTag_df.rename(columns=cols_to_rename, inplace=True)
        else:
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

        return attributeTag_df

    @write_df_excel
    def query_8(self) -> pd.DataFrame:
        # TODO:remove all tables except table_5_6_7_unioned
        # 1. get razor_protein_associated_with_peptide
        #         table_1,table_2,table_3,table_4,
        table_5_6_7_unioned = self.parse_MSGFjobs_MASIC_resultant()
        # 2. get annotations
        AnnotationSplitout_df = self.query_0()

        merged = table_5_6_7_unioned.merge(
            AnnotationSplitout_df,
            how="left",
            left_on=["RazorProtein"],
            right_on=["Protein"],
        )
        col_of_interest = [
            "PeptideSequence",
            "RazorProtein",
            "Product",
            "ec_number",
            "pfam",
            "ko",
            "cog",
        ]

        cols_to_rename = {"ec_number": "EC_Number", "ko": "KO", "cog": "COG"}
        PeptideRazorProteinAnnotated = merged[col_of_interest]
        PeptideRazorProteinAnnotated.rename(columns=cols_to_rename, inplace=True)

        return PeptideRazorProteinAnnotated

    @write_df_excel
    def get_FwdPeptideSequences(self) -> pd.DataFrame:
        FiltPeps_df = self.get_FiltPeps(self.resultant_df)
        FiltPeps_df["GeneCount"] = FiltPeps_df.groupby(["PeptideSequence"])[
            "Protein"
        ].transform("nunique")
        return FiltPeps_df[["PeptideSequence", "GeneCount"]]

    @write_df_excel
    def query_9(self) -> pd.DataFrame:
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

        return (
            merged[["PeptideSequence", "GeneCount", "Protein"]]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"])
        )

    def build_annotation_str(self, Protein: str, Product: str, pfam: str, ko: str, ec_number: str, cog: str) -> str:

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

    @write_df_excel
    def query_10(self) -> pd.DataFrame:
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

        return merged[
            ["PeptideSequence", "GeneCount", "Protein", "annotation"]
        ].sort_values(by=["PeptideSequence"])

    @write_df_excel
    def query_11(self) -> pd.DataFrame:
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

        return (
            table_10[["PeptideSequence", "GeneCount", "FullGeneList", "AnnotationList"]]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"])
        )

    @write_df_excel
    def query_12(self) -> pd.DataFrame:

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

        return (
            FiltPassPeptideAbundanceData[
                ["PeptideSequence", "SpectralCount", "sum(StatMomentsArea)"]
            ]
            .drop_duplicates()
            .sort_values(by=["PeptideSequence"])
        )

    @write_df_excel
    def query_13(self) -> pd.DataFrame:
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

        PeptideRazorProteinAnnotated_df = self.query_8()
        # saved for multiple call.
        self.query_8_result = PeptideRazorProteinAnnotated_df.copy()

        # left outer join 'PeptideSequence'
        merged_2 = merged_1.merge(
            PeptideRazorProteinAnnotated_df, how="left", on=["PeptideSequence"]
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

        return (
            merged_3[
                [
                    "DatasetName",
                    "PeptideSequence",
                    "RazorProtein",
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
            .drop_duplicates() # TODO check if this needs to be removed, probably does
            .sort_values(by=["DatasetName"])
        )

    @write_df_excel
    def query_14(self) -> pd.DataFrame:
        Peptide_Report = self.query_13()

        Peptide_Report["UniquePeptideCount"] = Peptide_Report.groupby(["RazorProtein"])[
            "PeptideSequence"
        ].transform("count")
        Peptide_Report["SummedSpectraCounts"] = Peptide_Report.groupby(["RazorProtein"])[
            "SpectralCount"
        ].transform("sum")
        Peptide_Report["SummedPeptideMASICAbundances"] = Peptide_Report.groupby(
            ["RazorProtein"]
        )["sum(MASICAbundance)"].transform("sum")

        return (
            Peptide_Report[
                [
                    "RazorProtein",
                    "UniquePeptideCount",
                    "SummedSpectraCounts",
                    "SummedPeptideMASICAbundances",
                ]
            ]
            .drop_duplicates() # TODO check if this needs to be removed, probably does
            .sort_values(by=["RazorProtein"])
        )

    @write_df_excel
    def query_15(self) -> pd.DataFrame:

        PeptideRazorProteinAnnotated_df = self.query_8_result  # self.query_8()
        temp_df = self.resultant_df.copy()
        temp_df["PeptideSequence"] = temp_df["Peptide"].str.extract(r"\.(.*)\.")
        merged = PeptideRazorProteinAnnotated_df.merge(
            temp_df,
            how="inner",
            left_on=["PeptideSequence"],
            right_on=["PeptideSequence"],
        )
        non_decoy_filtered_df = merged[~merged["Protein"].str.startswith("XXX")]

        return (
            non_decoy_filtered_df[["RazorProtein", "Protein"]]
            .drop_duplicates(subset=["RazorProtein", "Protein"]) # TODO check if this needs to be removed, probably does
            .sort_values(by=["RazorProtein"])
        )

    @write_df_excel
    def query_16(self) -> pd.DataFrame:
        RazorProteinGeneLists = self.query_15()  # .drop_duplicates(subset=['Protein'])

        # prepare protein list
        RazorProteinGeneLists["FullGeneList"] = RazorProteinGeneLists.groupby(
            ["RazorProtein"]
        )["Protein"].transform(lambda x: ", ".join(x))
        # sort them!
        get_list_func = lambda x: x.split(", ") if len(x.split(", ")) > 1 else x
        f = (
            lambda x: ", ".join(sorted(get_list_func(x)))
            if isinstance(get_list_func(x), list)
            else get_list_func(x)
        )
        RazorProteinGeneLists["FullGeneList"] = RazorProteinGeneLists[
            "FullGeneList"
        ].apply(f)

        RazorProteinGeneLists["GeneCount"] = RazorProteinGeneLists.groupby(
            ["RazorProtein"]
        )["Protein"].transform("count")
        RazorProteinGeneLists.sort_values(["FullGeneList"], ascending=True, inplace=True)

        return (
            RazorProteinGeneLists[["RazorProtein", "GeneCount", "FullGeneList"]]
            .drop_duplicates()
            .sort_values(by=["RazorProtein"])
        )

    def mod_build_annotation_str(
        self, RazorProtein, Protein, Product, pfam, ko, ec_number, cog
    ):

        if Product != Product:  # pd.isna(Product):
            # handle <class 'float'> nan
            return RazorProtein
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

    @write_df_excel
    def query_17(self) -> pd.DataFrame:

        AnnotationSplitout_df = self.query_0()

        table_15 = self.query_15()  # have RazorProtein
        ProteinRazorProteinsAnnotated = table_15.merge(
            AnnotationSplitout_df, how="left", left_on=["Protein"], right_on=["Protein"]
        )
        #         display(ProteinRazorProteinsAnnotated)
        ProteinRazorProteinsAnnotated["annotation"] = ProteinRazorProteinsAnnotated.apply(
            lambda x: self.mod_build_annotation_str(
                x.RazorProtein, x.Protein, x.Product, x.pfam, x.ko, x.ec_number, x.cog
            ),
            axis=1,
        )

        return (
            ProteinRazorProteinsAnnotated[["RazorProtein", "Protein", "annotation"]]
            .drop_duplicates()
            .sort_values(by=["RazorProtein"])
        )

    @write_df_excel
    def query_18(self) -> pd.DataFrame:

        RazorProteinAnnotationLists = self.query_17()

        RazorProteinAnnotationLists[
            "AnnotationList"
        ] = RazorProteinAnnotationLists.groupby(["RazorProtein"])["annotation"].transform(
            lambda x: " | ".join(x)
        )

        return (
            RazorProteinAnnotationLists[["RazorProtein", "AnnotationList"]]
            .drop_duplicates()
            .sort_values(by=["RazorProtein"])
        )

    @write_df_excel
    def query_19(self) -> pd.DataFrame:
        Peptide_Report_df = self.peptide_report.copy()[
            ["DatasetName", "RazorProtein", "Product", "EC_Number", "pfam", "KO", "COG"]
        ]
        RazorProteinGeneLists = self.query_16()
        merged_1 = Peptide_Report_df.merge(
            RazorProteinGeneLists,
            how="inner",
            left_on=["RazorProtein"],
            right_on=["RazorProtein"],
        )

        RazorProteinAnnotationLists = self.query_18()
        merged_2 = merged_1.merge(
            RazorProteinAnnotationLists,
            how="inner",
            left_on=["RazorProtein"],
            right_on=["RazorProtein"],
        )

        PeptideAbundanceInfo_df = self.query_14()
        merged_3 = merged_2.merge(
            PeptideAbundanceInfo_df,
            how="inner",
            left_on=["RazorProtein"],
            right_on=["RazorProtein"],
        )

        distinct_colm = [
            "DatasetName",
            "RazorProtein",
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

        return (
            merged_3[distinct_colm]
            .drop_duplicates(subset=distinct_colm)
            .sort_values(by=["DatasetName"])
        )


    @write_df_excel
    def query_20(self) -> pd.DataFrame:

        total_proteins_per_peptide = pd.DataFrame(
            {"total_protein_count": self.query_0()["Protein"].nunique()}, index=[0]
        )

        return total_proteins_per_peptide

    @write_df_excel
    def query_21(self) -> pd.DataFrame:
        temp_df = self.resultant_df.copy()
        filtered = self.get_FiltPeps_gen(temp_df)[["Scan", "Charge"]]
        PSM_Count = pd.DataFrame(
            {"PSM_Count": len(filtered.groupby(["Scan", "Charge"])["Scan"])}, index=[0]
        )

        return PSM_Count

    @write_df_excel
    def query_22(self) -> pd.DataFrame:
        temp_df = self.resultant_df.copy()
        Total_PSM_Count = pd.DataFrame(
            {
                "Total_PSM_Count": len(
                    temp_df[["Scan", "Charge"]].groupby(["Scan", "Charge"])["Scan"]
                )
            },
            index=[0],
        )

        return Total_PSM_Count

    @write_df_excel
    def query_23(self) -> pd.DataFrame:
        Table12 = pd.DataFrame(
            {
                "PSM_identification_rate": (
                    self.query_21().PSM_Count / self.query_22().Total_PSM_Count
                )
                * 100
            },
            index=[0],
        )

        return Table12

    @write_df_excel
    def query_24(self) -> pd.DataFrame:
        Peptide_Report_df = self.peptide_report.copy()
        unique_peptide_seq_count = (
            Peptide_Report_df.PeptideSequence.drop_duplicates().count()
        )
        RazorProtein_count = Peptide_Report_df.RazorProtein.drop_duplicates().count()
        mean_peptide_count = (unique_peptide_seq_count / RazorProtein_count)

        qc_metrics = pd.DataFrame(
            {
                "PSM_count": self.query_21().PSM_Count,
                "PSM_identification_rate": self.query_23().PSM_identification_rate,
                "unique_peptide_seq_count": unique_peptide_seq_count,
                "RazorProtein_count": RazorProtein_count,
                "mean_peptide_count": mean_peptide_count,
                "total_protein_count": self.query_20().total_protein_count
            },
            index=[0],
        )

        return qc_metrics

    @write_df_excel
    def create_peptide_report(self) -> pd.DataFrame:
        Peptide_Report_df = self.query_13()
        Peptide_Report_df.sort_values(["GeneCount"], ascending=False, inplace=True)

        return Peptide_Report_df

    @write_df_excel
    def create_protein_report(self) -> pd.DataFrame:
        Protein_Report_df = self.query_19()
        Protein_Report_df.sort_values(["GeneCount"], ascending=False, inplace=True)
        
        return Protein_Report_df

    @write_df_excel
    def create_qc_metrics(self) -> pd.DataFrame:
        qc_metrics_df = self.query_24()
        
        return qc_metrics_df

    def gen_reports(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        print(f"ReportGen start {self.dataset_id} : {self.faa_id}")
        self.peptide_report = self.create_peptide_report()
        protein_report = self.create_protein_report()
        qc_metrics_report = self.create_qc_metrics()
        print(f"ReportGen end {self.dataset_id} : {self.faa_id}")

        # rename to adhere nmdc.schema.json
        cols_to_rename = {
            "PeptideSequence": "peptide_sequence",
            "RazorProtein": "razor_protein",
            "FullGeneList": "all_proteins",
            "min(QValue)": "min_q_value",
            "SpectralCount": "peptide_spectral_count",
            "sum(MASICAbundance)": "peptide_sum_masic_abundance",
        }
        self.peptide_report.rename(columns=cols_to_rename, inplace=True)

        # rename to adhere nmdc.schema.json
        cols_to_rename = {
            "RazorProtein": "razor_protein",
            "FullGeneList": "all_proteins",
            "sum(MASICAbundance)": "peptide_sum_masic_abundance",
        }
        protein_report.rename(columns=cols_to_rename, inplace=True)

        # TODO: Change string type to "min_q_value", 'peptide_spectral_count', 'peptide_sum_masic_abundance' int!
        return self.peptide_report, protein_report, qc_metrics_report
    
    @staticmethod
    def parse_attributes(raw_attributes: str):
        if str(raw_attributes) == "nan":
            return None

        attributes = raw_attributes.split(";")
        attributes_map = {}
        for attr in attributes:
            if "=" in attr:
                key, value = attr.split("=", 1)
                attributes_map[key] = value
        return attributes_map

if __name__ == "__main__":

    fasta_txt_file = sys.argv[1]
    gff_file = sys.argv[2]
    resultant_file = sys.argv[3]
    dataset_id= sys.argv[4]
    faa_id = sys.argv[5]
    threshold= sys.argv[6]
    dataset_name= sys.argv[7]
    is_split_analysis= sys.argv[8]
    is_metagenome_free_analysis = sys.argv[9]

    print(f"{fasta_txt_file}\n{gff_file}\n{resultant_file}\n{dataset_id}\n{faa_id}\n{threshold}\n")

    is_split_analysis = is_split_analysis.rstrip().lower() == "true"
    is_metagenome_free_analysis = is_metagenome_free_analysis.rstrip().lower() == "true"

    data_obj = DataOutputtable(
        gff_file,
        resultant_file,
        fasta_txt_file,
        threshold,
        dataset_id,
        faa_id,
        dataset_name,
        is_split_analysis,
        is_metagenome_free_analysis
    )

    # data_obj.parse_MSGFjobs_MASIC_resultant()

    (
        peptide_report,
        protein_report,
        qc_metrics_report,
    ) = data_obj.gen_reports()

    # write dfs to file.
    peptide_report.to_csv(f"{dataset_id}_{faa_id}_Peptide_Report.tsv", sep="\t", index=False)
    protein_report.to_csv(f"{dataset_id}_{faa_id}_Protein_Report.tsv", sep="\t", index=False)
    qc_metrics_report.to_csv( f"{dataset_id}_{faa_id}_QC_metrics.tsv", sep="\t", index=False)
    
    # # flush and close excel writer
    # data_obj.writer.close()
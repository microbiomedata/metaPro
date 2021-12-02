def contaminate_fasta(self, files):

    contaminated_file = os.path.join(save_at, os.path.basename(files[0]))
    with open(contaminated_file, "w") as outfile:
        for file in files:
            with open(file) as infile:
                for line in infile:
                    outfile.write(line)
    return contaminated_file


if __name__ == "__main__":

    fasta_file = sys.argv[1]
    contaminant_file = sys.argv[2]


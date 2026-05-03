class RNASequence:
    CODON_TABLE = {
        "AUG": "Met",

        "UUU": "Phe", "UUC": "Phe",
        "UUA": "Leu", "UUG": "Leu",

        "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",

        "UAU": "Tyr", "UAC": "Tyr",
        "UAA": "STOP", "UAG": "STOP",

        "UGU": "Cys", "UGC": "Cys",
        "UGA": "STOP",
        "UGG": "Trp",

        "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",

        "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",

        "CAU": "His", "CAC": "His",
        "CAA": "Gln", "CAG": "Gln",

        "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",

        "AUU": "Ile", "AUC": "Ile", "AUA": "Ile",

        "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",

        "AAU": "Asn", "AAC": "Asn",
        "AAA": "Lys", "AAG": "Lys",

        "AGU": "Ser", "AGC": "Ser",
        "AGA": "Arg", "AGG": "Arg",

        "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",

        "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",

        "GAU": "Asp", "GAC": "Asp",
        "GAA": "Glu", "GAG": "Glu",

        "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
    }

    def __init__(self, sequence):
        self.__sequence = sequence
        self.validate_seq()

    def validate_seq(self):
        seq = self.__sequence.upper().replace(" ", "").replace("\n", "")
        valid_bases = {"A", "U", "G", "C"}

        if len(seq) == 0:
            raise ValueError("Sequence can't be empty.")

        for base in seq:
            if base not in valid_bases:
                raise ValueError(f"{base} not in {valid_bases}.")

        self.__sequence = seq

    def get_codons(self):
        seq = self.__sequence

        return [
            seq[i:i + 3]
            for i in range(0, len(seq), 3)
            if len(seq[i:i + 3]) == 3
        ]

    def translate(self):
        protein = []

        start = self.get_start_codon()

        if start == -1:

            return protein

        for i in range(start, len(self.__sequence), 3):

            codon = self.__sequence[i:i + 3]

            if len(codon) < 3:

                break

            amino = RNASequence.CODON_TABLE.get(codon)

            if amino == "STOP":

                break

            if amino is not None:

                protein.append(amino)

        return protein
    def get_start_codon(self):
        return self.__sequence.find("AUG")

    def has_start_codon(self):
        return self.get_start_codon() != -1

    def get_stop_codon(self):
        seq = self.__sequence
        stop_codons = {"UAA", "UAG", "UGA"}

        start_codon = self.get_start_codon()

        if start_codon == -1:
            return -1

        for i in range(start_codon + 3, len(seq), 3):
            codon = seq[i:i + 3]

            if len(codon) < 3:
                break

            if codon in stop_codons:
                return i

        return -1

    def has_stop_codon(self):
        return self.get_stop_codon() != -1

    def get_coding_region(self):
        start = self.get_start_codon()
        stop = self.get_stop_codon()

        if start == -1 or stop == -1:
            return ""

        return self.__sequence[start:stop]

    def amino_acids_count(self):
        protein = self.translate()
        amino_dict = {}

        for amino in protein:
            if amino not in amino_dict:
                amino_dict[amino] = 1
            else:
                amino_dict[amino] += 1

        return amino_dict

    def predict_protein_functionality(self):
        protein = self.translate()

        if not self.has_start_codon():
            return "Start codon not found. RNA sequence is unlikely to translate properly."

        if not self.has_stop_codon():
            return "Stop codon not found. Protein may be incomplete or extended."

        if len(protein) == 0:
            return "No amino acid product detected."

        if len(protein) < 30:
            return "Protein is less than 30 amino acids. It may be a small peptide or nonfunctional fragment."

        return (
            "Sequence contains a start codon, stop codon, and meets basic length requirements.\n"
            "It will likely produce a functional protein."
        )

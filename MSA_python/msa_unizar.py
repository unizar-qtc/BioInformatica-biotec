#!/usr/bin/env python3

"""
=======================================================================
  MSA - UNIZAR
=======================================================================

  Python3 module with Multiple Sequence Alignment functions

  Bioinformatics Laboratory
  Biotechnology Degree @ University of Zaragoza (Spain)


  Functions
  ---------
    download_pdb(pdb_id)
    read_sequence(pdb_id, pdb_chain)
    get_alignment_file(pdb_id, pdb_chain, aa_seq)
    consensus_alignment(pdb_id, pdb_chain, aa_seq, muscle_exe)

"""


##  DEPENDENCIES  #####################################################

from io import StringIO
from urllib import request

import Bio
from Bio import AlignIO, SeqIO, pairwise2
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62

try:
    import colorama
    from colorama import Fore, Style
except ImportError:
    colorama = None


##  FUNCTIONS  ########################################################

def download_pdb(pdb_id:str) -> None:
    """
        Download a PDB file from rcsb.org based on the 4 letter/numbers ID

        Parameters
        ----------
        pdb_id : str
            PDB ID (4 letters/numbers, case sensitive)
    """

    file_name = pdb_id + ".pdb"
    request.urlretrieve('https://files.rcsb.org/download/' + file_name, file_name)

def read_sequence(pdb_id:str, pdb_chain:str) -> 'Bio.Seq.Seq':
    """
        Read the annotated sequence of a chain in a PDB file

        Parameters
        ----------
        pdb_id : str
            PDB ID (4 letters/numbers, case sensitive)
        pdb_chain : str
            chain ID (one letter)

        Returns
        -------
        aminoacid_seq : Bio.Seq.Seq
            amino acid sequence
    """

    structure_file = pdb_id + ".pdb"
    pdb_chain      = pdb_chain.upper()

    myseq_full = SeqIO.parse(structure_file, "pdb-seqres")
    for chain in myseq_full:
        if chain.id[5:] == pdb_chain:
            aminoacid_seq = chain.seq   # Extract the amino acid sequence of the protein chain of interest
            return aminoacid_seq

def get_alignment_file(pdb_id:str, pdb_chain:str, aa_seq:'Bio.Seq.Seq') -> None:
    """
        Obtain a BLAST alignment from an amino acid sequence.

        Several files are generated and named based on the PDB ID and
        chain ID (e.g. XXXX, Y):
            blast_XXXX.xml
            blast_XXXX_clean.xml
            XXXX-Y_in.fasta
            XXXX-Y_align.fasta

        Parameters
        ----------
        pdb_id : str
            PDB ID (4 letters/numbers, case sensitive)
        pdb_chain : str
            chain ID (one letter)
        aa_seq : Bio.Seq.Seq
            amino acid sequence
    """

    file_out_name = f"blast_{pdb_id}.xml"
    file_out = open(file_out_name, "w")
    print("BLASTing....")
    handle = NCBIWWW.qblast("blastp", "nr", aa_seq, hitlist_size=10000, gapcosts="10 1", expect=0.0001)     # With nonredundant database
    file_out.write(handle.read())
    file_out.close()
    _clean_xml(file_out_name)

    handle = open(f"{file_out_name[:-4]}_clean.xml")
    blast_record = NCBIXML.read(handle)                                 # (handle) is used as an argument without saving a BLAST file
    ids = ["my_protein"]                                                # The given sequence is the only one without ID in BLAST, so this one is given
    sequences = [str(aa_seq)]                                           # Further sequence clustering and alignments need to be performed
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if (hsp.identities/float(hsp.align_length)) > 0.7 and "Z" not in hsp.sbjct and "B" not in hsp.sbjct:
                this_seq = ""
                for char in hsp.sbjct:
                    if char != "-":
                        this_seq += char
                sequences.append(this_seq)                              # Takes the other sequences (the aligned part) and gets it into a list
                name = alignment.title
                while name[0] == "|":
                    name = name[1:]
                if name in ids:
                    name = "A" + name
                    if name == "AAAAAAAAAAA":
                        break
                ids.append(name)                                        # Takes the title of the other sequences
    short_ids = []                                                      # Empty lists for the records' properties
    desc_list = []
    new_sequences = []
    count = 0
    for item in ids:
        desc = "NA"                                                     # For proteins without description, NA will be theirs (only expected in "my protein")
        if "|" in item:
            index = 0
            for i in range(2):
                index = item.find("|", index + 1, len(item))            # For getting the short id (just until the fourth bar)
            desc = item[index + 1:]                                     # The rest after (+1) the fourth bar is the description
            item = item[3:index]                                        # The three first characters indicate the source, if not deleted, it may lead to
                                                                        # errors in writing strict Phylip such as getting two identical IDs from two
                                                                        # different sequences.
            j = 0
            while len(item) <= 10:
                if j < len(desc):
                    item += desc[j]
                    j += 1
                else:
                    break
        if "XP_" in item or "NP_" in item or "WP_" in item or "YP_" in item:
            index_n = item.find("NP_")
            index_x = item.find("XP_")
            index_w = item.find("WP_")
            index_y = item.find("YP_")
            index_1 = -5
            index_list = [index_n, index_x, index_w, index_y]
            while index_1 < 0:
                index_1 = min(index_list)
                index_list[index_list.index(index_1)] = 999
            if "." in item:
                index_2 = 0
                while index_2 < index_1:
                    index_2 = item.find(".")
                else:
                    item = item[index_1 + 3:index_2] + item[index_2 + 1:]
            else:
                item = item[index_1 + 3:]
        while item[0] == "|":
            item = item[1:]
        if item not in short_ids:
            short_ids.append(item)
            desc_list.append(desc)
            new_sequences.append(sequences[count])
        count += 1
    fasta_list = []                                                     # To keep all sequences in a list.
    for i in range(len(new_sequences)):
        my_seq = SeqRecord(Seq(new_sequences[i]))                       # Make a SeqRecord with each sequence
        my_seq.id = short_ids[i]                                        # Add their properties to each SeqRecord (they are the same order in the list)
        my_seq.description = desc_list[i]
        fasta_list.append(my_seq)                                       # Append SeqRecord objects to a new list, for writing them
    doc_handle = f"{pdb_id}-{pdb_chain}_in.fasta"
    align_handle = f"{pdb_id}-{pdb_chain}_align.fasta"
    SeqIO.write(fasta_list, doc_handle, "fasta")                        # Writes the SeqRecord objects as sequences in the file doc_handle
    if len(fasta_list) >= 250:
        SeqIO.write(fasta_list[0:250], align_handle, "fasta")
    else:
        SeqIO.write(fasta_list[0:len(fasta_list)], align_handle, "fasta")

def consensus_alignment(pdb_id:str, pdb_chain:str, aa_seq:'Bio.Seq.Seq', muscle_exe:str) -> None:
    """
        Get a consensus sequence based on a BLAST result and print on
        the screen its alignment to the original sequence

        Parameters
        ----------
        pdb_id : str
            PDB ID (4 letters/numbers, case sensitive)
        pdb_chain : str
            chain ID (one letter)
        aa_seq : Bio.Seq.Seq
            amino acid sequence
        muscle_exe : str
            MUSCLE executable (relative path)
    """

    file_in = f"{pdb_id}-{pdb_chain}_align.fasta"
    muscle_cline = MuscleCommandline(muscle_exe, input=file_in)
    stdout, stderr = muscle_cline()
    alignment = AlignIO.read(StringIO(stdout), "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)            # Object for studying properties of the alignment
    consensus = summary_align.dumb_consensus(threshold=0.5)     # Makes the simple consensus, with X as the no-consensus symbol
    res = pairwise2.align.globalds(aa_seq, consensus, blosum62, -10, -0.5)
    seq1 = res[0][0]
    seq2 = res[0][1]
    _fancy_seq_print(seq1, seq2)


##  PRIVATE FUNCTIONS  ################################################

def _clean_xml(file_dir:str) -> None:
    """
        Cleans a XML file, creating a new one starting with the same name
        and ending as '_clean.xml'

        Parameters
        ----------
        file_dir : str
            xml file (relative or absolute path)
    """

    file_in  = open(file_dir)
    file_out = open(".".join(file_dir.split(".")[:-1])+"_clean.xml", "w")

    i = 99
    check_dict = {"IH": False, "IS": False, "It": False, "BOIt": False, "BO": False}
    for line in file_in.readlines():
        if "CREATE_VIEW" in line:
            i = 0
        elif "</Iteration_hits>" in line:
            check_dict["IH"] = True
            file_out.write(line)
        elif "<Iteration_stat>" in line:
            check_dict["IS"] = True
            file_out.write(line)
        elif "</Iteration>" in line:
            check_dict["It"] = True
            file_out.write(line)
        elif "</BlastOutput_iterations>" in line:
            check_dict["BOIt"] = True
            file_out.write(line)
        elif "</BlastOutput>" in line:
            check_dict["BO"] = True
            file_out.write(line)
        elif i > 2:
            file_out.write(line)
        i += 1
    file_in.close()
    if not check_dict["IH"]:
        file_out.write("</Iteration_hits>\n")
    if not check_dict["IS"]:
        file_out.write("  <Iteration_stat>\n"
                       "    <Statistics>\n"
                       "      <Statistics_db-num>302672228</Statistics_db-num>\n"
                       "      <Statistics_db-len>1901495170</Statistics_db-len>\n"
                       "      <Statistics_hsp-len>0</Statistics_hsp-len>\n"
                       "      <Statistics_eff-space>0</Statistics_eff-space>\n"
                       "      <Statistics_kappa>0.024</Statistics_kappa>\n"
                       "      <Statistics_lambda>0.243</Statistics_lambda>\n"
                       "      <Statistics_entropy>0.1</Statistics_entropy>\n"
                       "    </Statistics>\n"
                       "  </Iteration_stat>\n")
    if not check_dict["It"]:
        file_out.write("</Iteration>\n")
    if not check_dict["BOIt"]:
        file_out.write("</BlastOutput_iterations>\n")
    if not check_dict["BO"]:
        file_out.write("</BlastOutput>\n")
    file_out.close()

def _fancy_seq_print(seq1:'Bio.Seq.Seq', seq2:'Bio.Seq.Seq', line_lenght:int=50) -> None:
    """
        Print a sequence alignment and comparison in a fancy way

        Parameters
        ----------
        seq1 : Bio.Seq.Seq
            amino acid sequence of reference
        seq2 : Bio.Seq.Seq
            amino acid sequence for comparison
        line_length : int, optional
            number of residues to print per line
    """

    if colorama:
        colorama.init()
        colors = [Fore.RED if j=='X' else Fore.CYAN if i==j else Fore.YELLOW for i, j in zip(seq1, seq2)]
        c_normal = Style.RESET_ALL
    else:
        colors = [""]*len(seq1)
        c_normal = ""

    for i in range(len(seq1)):
        if i % line_lenght == 0:
            print(f"\n\n{i+1:4d}  {seq1[i:min(len(seq1),i+line_lenght)]}\n"+" "*6, end="")
        print(f"{colors[i]}{seq2[i]}{c_normal}", end="")
    print()

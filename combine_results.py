import os
import sys
import shutil
from Bio import SeqIO

def combine_results(msa_file, db_files, out_file):
    # msa_file: alignment file path
    # db_files: a list of structure files in dot-bracket format
     
    # read structures
    db_struc = {}
    for file in db_files:
        file_name = os.path.basename(file)
        idx = int(file_name.split("_")[0])
        with open(file, "r") as f:
            line = f.readlines()[-1].strip()
        db_struc[idx - 1] = line

    # read alignment
    msa_data = {}
    names = []
    for record in SeqIO.parse(msa_file, "fasta"):
        name = record.id
        names.append(name)
        msa_data[name] = str(record.seq)

    # align structures
    align_strucs = {}
    for i, name in enumerate(names):
        struc = db_struc[i]
        align_struc = ""
        align_seq = msa_data[name]

        pos = 0
        for nuc in align_seq:
            if nuc == "-":
                align_struc += "-"
            else:
                align_struc += struc[pos]
                pos += 1

        assert len(align_struc) == len(align_seq)
        align_strucs[name] = align_struc

    # output 
    with open(out_file, "w") as f:
        for name in names:
            f.write(">%s\n" % name)
            f.write("%s\n" % msa_data[name])
            f.write("%s\n" % align_strucs[name])


if __name__ == "__main__":
    out_dir = sys.argv[1]
    out_file = "ltf.out"

    # catch error info when combining results
    try:
        files = os.listdir(out_dir)
        db_files = [os.path.join(out_dir, file) for file in files if file.endswith(".db")]
        msa_file = "%s/output.aln" % out_dir

        combine_results(msa_file, db_files, out_file)
    except OSError as e:
        print("Error: %s : %s" % ("combine results", e.strerror))

    # optional: remove linearturbofold output dir
    # try:
    #     shutil.rmtree(out_dir)
    # except OSError as e:
    #     print("Error: delete %s : %s" % (out_dir, e.strerror))

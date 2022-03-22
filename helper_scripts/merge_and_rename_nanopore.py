from glob import glob
import dinopy
import sys


def merge_and_rename(folder_path):
    barcode = folder_path.split("/")[-2]
    print("barcode: " + barcode)
    with dinopy.FastqWriter(folder_path + "/" + barcode + "_A_R1.fastq.gz") as out:
        for path in glob(folder_path + "/*"):
            inp = dinopy.FastqReader(path)
            for entry in inp.reads(quality_values=True):
                name = entry.name.decode().split(" ")[0].split("-")[0] + "_" + barcode[7:]
                out.write(entry.sequence, name.encode(), entry.quality)


if __name__ == "__main__":
        f_path_root = sys.argv[1]
        for f_path in glob(f_path_root):
            merge_and_rename(f_path)

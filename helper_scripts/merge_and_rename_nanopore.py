from glob import glob
import dinopy
import sys
import os


def merge_and_rename(folder_path):
    barcode = folder_path.split("/")[-2]
    print("barcode: " + barcode)
    if not os.path.exists("/".join(folder_path.split("/")[:-3]) + "/" + barcode + "_A_R1.fastq.gz"):
        with dinopy.FastqWriter("/".join(folder_path.split("/")[:-3]) + "/" + barcode + "_A_R1.fastq.gz") as out:
            for path in glob("/".join(folder_path.split("/")[:-1]) + "/*"):
                inp = dinopy.FastqReader(path)
                for entry in inp.reads(quality_values=True):
                    name = entry.name.decode().split(" ")[0].split("-")[0] + "_" + barcode[7:]
                    out.write(entry.sequence, name.encode(), entry.quality)


if __name__ == "__main__":
        f_path_root = sys.argv[1] # path/to/barcodefolders/*/*.fastq.qz
        for f_path in glob(f_path_root):
            merge_and_rename(f_path)

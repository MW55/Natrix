import sys
import yaml
import glob
import os

yaml_file = sys.argv[1]
with open(yaml_file) as f_:
    config = yaml.load(f_, Loader=yaml.FullLoader)

fileslist = [os.path.basename(x) for x in glob.glob(config['general']['filename'] + '/*.fastq.gz')]

if len(fileslist) == 0:
    print("No Files found to symlink..exiting program.")
    sys.exit()

#getting an example filename to process
sample = fileslist[0]

#remove fastq.gz
suffix = '.fastq.gz'
sample = sample[:-len(suffix)]

print(f"Using \'{sample}\' as an example to define how to read the provided fastq.gz files.")

#identify separator used
separator = None
while not separator:
    separator = input("Please input the separator used:")

splitted_sample = sample.split(separator)
print("The following file parts have been indentified:")
for (i, val) in enumerate(splitted_sample):
    print(i, val)

#identify sample name parts
name_index = int(input("Please input the number corresponding to the sample name:"))

#identify unit part if split sample
unit_index= -1
if config['merge']['filter_method'] == 'split_sample':
    unit_index = int(input("Please input the number corresponding to the unit sample:"))

#identify read part if paired_end
read_index= -1
if config['merge']['paired_End']:
    read_index = int(input("Please input the number corresponding to the forward or reverse read of the sample:"))

#mapping of files to roughly natrix schema
def get_mapped_filename(file):
    splitted = file[:-len(suffix)].split(separator)
    mapped_file = splitted[name_index]
    mapped_file += '_'
    if config['merge']['filter_method'] == 'split_sample':
        mapped_file += splitted[unit_index]
    else:
        mapped_file += 'A'
    mapped_file += '_'
    if config['merge']['paired_End']:
        mapped_file += splitted[read_index]
    else:
        mapped_file += 'R1'
    mapped_file += suffix
    return mapped_file

mapped_files = list(map(get_mapped_filename, fileslist))

#remapping units
unit_list = [x[:-len(suffix)].split('_')[1] for x in mapped_files]
unit_list = list(set(unit_list))
if(len(unit_list) > 2):
    print("Too many units per sample?..exiting program.")
    sys.exit()
if(len(unit_list) == 2):
    print("The following units have been indentified:")
    for (i, val) in enumerate(unit_list):
        print(i, val)
    A = int(input("Please input the number corresponding to the 'A' unit:"))
else:
    A = 0
def remap_units(file):
    splitted = file[:-len(suffix)].split('_')
    splitted[1] = 'A' if splitted[1] == unit_list[A] else 'B'
    return '_'.join(splitted)+suffix
mapped_files = list(map(remap_units, mapped_files))

#remapping reads
read_list = [x[:-len(suffix)].split('_')[2] for x in mapped_files]
read_list = list(set(read_list))
if(len(unit_list) > 2):
    #should this happen?
    print("Too many reads per sample?..exiting program.")
if(len(unit_list) == 2):
    print("The following reads have been indentified:")
    for (i, val) in enumerate(read_list):
        print(i, val)
    forward = int(input("Please input the number corresponding to the forward read:"))
else:
    forward = 0
def remap_reads(file):
    splitted = file[:-len(suffix)].split('_')
    splitted[2] = 'R1' if splitted[2] == read_list[forward] else 'R2'
    return '_'.join(splitted)+suffix
mapped_files = list(map(remap_reads, mapped_files))

#check mapping
filename_len = len(max(fileslist, key=len))
print("Files will be mapped as follows:")
for old_file,new_file in zip(fileslist,mapped_files):
    print(f"{old_file.ljust(filename_len)} => {new_file}")
if len(mapped_files) != len(set(mapped_files)):
    print("Some files would map to the same name..exiting program.")
    sys.exit()

print("Writing configuration..")


config['sample'] = dict()

config['sample'] = {
    'separator': separator,
    'name_idx': name_index,
    'unit_idx': unit_index,
    'unit_A_identifier': unit_list[A],
    'read_idx': read_index,
    'read_forward_identifier': read_list[forward],
}


with open(yaml_file, 'w') as outfile:
    yaml.dump(config, outfile, default_flow_style=False)

print('..done')

# # symlinking actual symlinking done in pipeline
# if not os.path.exists('symlink_test'):
#     os.mkdir('symlink_test')
# print("symlinking files to symlink_test/")
# for old_file,new_file in zip(fileslist,mapped_files):
#     os.symlink(config['general']['filename'] + '/'+old_file, 'symlink_test/'+new_file)
# print("symlinks done.")

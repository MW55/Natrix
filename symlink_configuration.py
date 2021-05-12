import sys
import yaml
import glob
import os

yaml_file = sys.argv[1]
with open(yaml_file) as f_:
    config = yaml.load(f_, Loader=yaml.FullLoader)

fileslist = [os.path.basename(x) for x in glob.glob(config['general']['filename'] + '/*.fastq.gz')]
suffix = '.fastq.gz'

if len(fileslist) == 0:
    print("No Files found to symlink. Exiting.")
    sys.exit()

if ('sample' in config):
    print(f"""The current symlink configuration is:
separator: {config['sample']['separator']}
name_idx: {config['sample']['name_idx']}
unit_idx: {config['sample']['unit_idx']}
read_idx: {config['sample']['read_idx']}
unit_A_identifier: {config['sample']['unit_A_identifier']}
read_forward_identifier: {config['sample']['read_forward_identifier']}
""")
    while True:
        new_conf = input("Do you want to change this configuration? (Enter yes or no):\n")
        if new_conf.lower() not in ['yes', 'no']:
            print('Invalid Input. Please answer with yes or no')
        else:
            break
else:
    new_conf = 'yes'

if new_conf == 'no':
    while True:
        check_conf = input(
            "Do you want to check the current configuration against the input files? (Enter yes or no):\n")
        if check_conf.lower() not in ['yes', 'no']:
            print('Invalid Input. Please answer with yes or no')
        else:
            break
    if check_conf == 'no':
        print("Nothing else to do. Exiting.")
        sys.exit()
    separator = config['sample']['separator']
    name_index = config['sample']['name_idx']
    unit_index = config['sample']['unit_idx']
    read_index = config['sample']['read_idx']
    A = config['sample']['unit_A_identifier']
    forward = config['sample']['read_forward_identifier']
else:
    # getting an example filename to process
    sample = fileslist[0]

    # remove fastq.gz
    sample = sample[:-len(suffix)]

    print(f"Using \'{sample}\' as an example to define how to read the provided sample files.")

    # identify separator used
    while True:
        separator = input("Please input the separator used:\n")
        if separator == "":
            print('Invalid Input.')
        else:
            break

    splitted_sample = sample.split(separator)
    print("The following file parts have been indentified:")
    for (i, val) in enumerate(splitted_sample):
        print(i, val)

    # identify sample name parts
    while True:
        try:
            name_index = int(input("Please input the number corresponding to the sample name:\n"))
        except ValueError:
            print("Invalid input. Please enter a valid number.")
            continue
        else:
            if name_index >= len(splitted_sample) or name_index < 0:
                print("Invalid input. Please enter a valid number.")
            else:
                break

    # identify unit part if split sample
    if config['merge']['filter_method'] == 'split_sample':
        print("Sample parts:")
        for (i, val) in enumerate(splitted_sample):
            print(i, val)
        while True:
            try:
                unit_index = int(input("Please input the number corresponding to the sample unit:\n"))
            except ValueError:
                print("Invalid input. Please enter a valid number.")
                continue
            else:
                if unit_index >= len(splitted_sample) or unit_index < 0:
                    print("Invalid input. Please enter a valid number.")
                else:
                    break
        unit_options = sorted(list(set(map(lambda x: x[:-len(suffix)].split(separator)[unit_index], fileslist))))
        print("The following units have been indentified:")
        for (i, val) in enumerate(unit_options):
            print(i, val)
        if len(unit_options) > 2:
            print("Too many different units for split-sample approach. Exiting.")
            sys.exit()
        if len(unit_options) < 2:
            print("Too few different units for split-sample approach. Exiting.")
            sys.exit()
        else:
            while True:
                try:
                    A_index = int(input("Please input the number corresponding to the 'A' unit:\n"))
                except ValueError:
                    print("Invalid input. Please enter a valid number.")
                    continue
                else:
                    if A_index >= len(unit_options) or A_index < 0:
                        print("Invalid input. Please enter a valid number.")
                    else:
                        A = unit_options[A_index]
                        break
    else:
        A = 'A'
        unit_index = -1

    # identify read part if paired_End
    if config['merge']['paired_End']:
        print("Sample parts:")
        for (i, val) in enumerate(splitted_sample):
            print(i, val)
        while True:
            try:
                read_index = int(
                    input("Please input the number corresponding to the forward or reverse read of the sample:\n"))
            except ValueError:
                print("Invalid input. Please enter a valid number.")
                continue
            else:
                if read_index >= len(splitted_sample) or read_index < 0:
                    print("Invalid input. Please enter a valid number.")
                else:
                    break
        read_options = sorted(list(set(map(lambda x: x[:-len(suffix)].split(separator)[read_index], fileslist))))
        print("The following units have been indentified:")
        for (i, val) in enumerate(read_options):
            print(i, val)
        if len(read_options) > 2:
            print("Too many different reads for paired end. Exiting.")
            sys.exit()
        if len(read_options) < 2:
            print("Too few different reads for paired end. Exiting.")
            sys.exit()
        else:
            while True:
                try:
                    forward_index = int(input("Please input the number corresponding to the forward read:\n"))
                except ValueError:
                    print("Invalid input. Please enter a valid number.")
                    continue
                else:
                    if forward_index >= len(read_options) or forward_index < 0:
                        print("Invalid input. Please enter a valid number.")
                    else:
                        forward = read_options[forward_index]
                        break
    else:
        forward = 'R1'
        read_index = -1


# check mapping
def map_file(file):
    splitted = file[:-len(suffix)].split(separator)
    mapped = splitted[name_index].replace('_', '-')
    mapped += '_'
    mapped += 'A' if config['merge']['filter_method'] == 'not_split' or A == splitted[unit_index] else 'B'
    mapped += '_'
    mapped += 'R1' if not config['merge']['paired_End'] or forward == splitted[read_index] else 'R2'
    mapped += suffix
    return mapped


mapped_files = {f: map_file(f) for f in fileslist}

filename_len = len(max(fileslist, key=len))
print("Files will be mapped as follows:")
for old_file, new_file in zip(fileslist, mapped_files):
    print(f"{old_file.ljust(filename_len)} => {new_file}")
if len(mapped_files) != len(set(mapped_files)):
    print("Some files would map to the same name..exiting program.")
    sys.exit()

if new_conf == 'yes':
    print("Saving configuration..")

    config['sample'] = dict()

    config['sample'] = {
        'separator': separator,
        'name_idx': name_index,
        'unit_idx': unit_index,
        'unit_A_identifier': A,
        'read_idx': read_index,
        'read_forward_identifier': forward,
    }

    with open(yaml_file, 'w') as outfile:
        yaml.dump(config, outfile, default_flow_style=False, sort_keys=False)

    print('..saving done.')
print('Nothing else to do. Exiting.')

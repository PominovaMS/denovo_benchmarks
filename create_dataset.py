import argparse
import os
from tqdm import tqdm
from oktoberfest.runner import run_job
from dataset_utils import *
from dataset_config import get_config


Q_VAL_THRESHOLD = 0.01
MAX_SPECTRA_PER_FILE = 20000

# Paths parsing
parser = argparse.ArgumentParser()
parser.add_argument("--config_path", type=str, help="path to the dataset config file")
args = parser.parse_args()

# Config setup
config = get_config(args.config_path)
dset_id = config.download.dset_id
dset_name = config.name

# Output dirs setup
for data_dir in [
    RAW_DATA_DIR, 
    MZML_DATA_DIR, 
    RESCORED_DATA_DIR, 
    DATASET_STORAGE_DIR
]:
    os.makedirs(data_dir, exist_ok=True)

# Create dirs for intermediate dataset files
raw_files_dir = os.path.join(RAW_DATA_DIR, dset_id)
mzml_files_dir = os.path.join(MZML_DATA_DIR, dset_name) # TODO: rename to smth search_files_dir? 
rescored_files_dir = os.path.join(RESCORED_DATA_DIR, dset_name)
mgf_files_dir = os.path.join(DATASET_STORAGE_DIR, dset_name, "mgf")

for data_dir in [
    raw_files_dir, mzml_files_dir, rescored_files_dir, mgf_files_dir,
]:
    os.makedirs(data_dir, exist_ok=True)

print("Creating dataset at:", mgf_files_dir)
fnames = [os.path.splitext(f)[0] for f in os.listdir(mgf_files_dir)]
print("Existing files:")
print(fnames)
# TODO: list summary of csv with labels instead?

# Get list of files to process 
files_list = get_files_list(config.download)
print("Processing files:\n", list(files_list.keys()))

# If already have pepxml, no need to run DB search (TODO: add flag for FORCED RE-RUN)
if not set(fname.lower() + ".pin" for fname in files_list).issubset(
        set(fname.lower() for fname in os.listdir(mzml_files_dir))
    ):
    # Prepare mzml (or other DB search input format) files if they don't exist
    if not set(fname.lower() + config.db_search.ext for fname in files_list).issubset(
            set(fname.lower() for fname in os.listdir(mzml_files_dir))
        ):
        # Download raw files from repository if needed
        if not all(
            os.path.exists(os.path.join(raw_files_dir, file_path)) for file_path in files_list.values()
        ):
            download_files(config.download, files_list)

        # Unpack if not unpacked yet
        if config.download.ext.endswith(".zip"):
            base_file_ext = config.download.ext[:-len(".zip")]
            # If files can be directly processed by MSFragger (.d, .mzml), unpack to mzml_files_dir
            if base_file_ext == ".d":
                unpack_dir = mzml_files_dir
            else:
                unpack_dir = raw_files_dir

            if not set(fname.lower() + base_file_ext for fname in files_list).issubset(
                set(fname.lower() for fname in os.listdir(unpack_dir))
            ):
                print(f"Files will be unpacked to {unpack_dir}")
                fnames = list(files_list.values())
                for fname in fnames:
                    file_path = os.path.join(raw_files_dir, fname)
                    shutil.unpack_archive(filename=file_path, extract_dir=unpack_dir)

        if config.download.ext in [".raw", ".wiff"]:
            # TODO: add: or (config.download.ext == ".d.zip" and config.search.ext == ".mzml")
            # Convert raw files to mzml (msconvert) # TODO: check for files from files_list only?
            convert_raw(dset_id, files_list, mzml_files_dir, target_ext=config.db_search.ext) # also always mzml??

        else:
            print(f"Unknown file extension {config.download.ext}")
    
    # Generate decoys for DB search
    db_w_decoys_path = generate_decoys_fasta(dset_name, config.db_search.database_path)

    # Run DB search (MSFragger)
    if config.db_search.n_db_splits == 1:
        run_database_search(dset_name, db_w_decoys_path, config.db_search)
    else:
        run_database_search_split(dset_name, db_w_decoys_path, config.db_search)

# Rescore DB search results (MSBooster + Percolator)
# Create features with MSBooster
# If already have _rescore.pin, no need to run feature prediction (TODO: add flag for FORCED RE-RUN)
file_prefix = "rescore"
# [! uncomment to use deep learning-based features in rescoring]
if not set(fname.lower() + f"_{file_prefix}.pin" for fname in files_list).issubset(
        set(fname.lower() for fname in os.listdir(mzml_files_dir))
    ):
    # TODO: ideally should also check whether mzml files exist and add them otherwise
    get_psm_rescoring_features(dset_name, config.rescoring)

# Run Percolator on created features
# TODO: any checks whether PSMs rescoring results already available?
if not f"{file_prefix}.percolator.psms.txt" in os.listdir(rescored_files_dir):
    run_psm_rescoring(dset_name, config.rescoring, files_list)

# Prepare unlabeled mgf files if they don't exist
if not set(fname.lower() + ".mgf" for fname in files_list).issubset(
        set(fname.lower() for fname in os.listdir(mgf_files_dir))
    ):
    # Download raw files from repository if needed
    if not all(os.path.exists(os.path.join(raw_files_dir, file_path)) for file_path in files_list.values()):
        download_files(config.download, files_list)

    # Unpack if not unpacked yet
    if config.download.ext.endswith(".zip"):
        base_file_ext = config.download.ext[:-len(".zip")]
        # If files can be directly processed by MSFragger (.d, .mzml), unpack to mzml_files_dir
        if base_file_ext == ".d":
            unpack_dir = mzml_files_dir
        else:
            unpack_dir = raw_files_dir

        if not set(fname.lower() + base_file_ext for fname in files_list).issubset(
            set(fname.lower() for fname in os.listdir(unpack_dir))
        ):
            print(f"Files will be unpacked to {unpack_dir}")
            fnames = list(files_list.values())
            for fname in fnames:
                file_path = os.path.join(raw_files_dir, fname)
                shutil.unpack_archive(filename=file_path, extract_dir=unpack_dir)

    if config.download.ext in [".raw", ".wiff"]:
        # Convert raw files to mgf (msconvert)
        convert_raw(dset_id, files_list, mgf_files_dir, target_ext=".mgf")

    elif config.download.ext == ".d.zip":
        # After DB search files should be already created 
        # from uncalibrated.mgf produced by MSFragger
        print("No mgf files found for .d data")
    
    else:
        print(f"Unknown file extension {config.download.ext}")


# Load DB search + rescoring results
results_path = os.path.join(rescored_files_dir, f"{file_prefix}.percolator.psms.txt")
results_df = pd.read_csv(results_path, sep="\t")
results_df = results_df[results_df["q-value"] < Q_VAL_THRESHOLD][["PSMId", "peptide", "q-value"]]
results_df["title"] = results_df["PSMId"].apply(lambda x: "_".join(x.split("_")[:-1]))

# filter and keep only spectra with defined charge
# get 0-based idxs from mgf files
spectra_idxs_0 = []

for fname in tqdm(files_list):
    input_path = os.path.join(mgf_files_dir, fname + ".mgf")
    spectra = mgf.read(input_path)
    
    spectra_filtered = [
        spectrum for spectrum in tqdm(spectra)
        if spectrum is not None and "charge" in spectrum["params"]
    ]

    # split filtered spectra into multiple files if needed
    n_splits = len(spectra_filtered) // MAX_SPECTRA_PER_FILE + int(
        len(spectra_filtered) % MAX_SPECTRA_PER_FILE > 0
    )
    print(f"Split {fname} file into {n_splits} files")
    
    for k in range(n_splits):
        start_idx = k * MAX_SPECTRA_PER_FILE
        end_idx = min((k + 1) * MAX_SPECTRA_PER_FILE, len(spectra_filtered))
        spectra_chunk = spectra_filtered[start_idx:end_idx]

        output_fname = f"{fname}_{k}"
        output_path = os.path.join(mgf_files_dir, output_fname + ".mgf")
        # TODO: add flag to force RE-WRITING
        if not os.path.exists(output_path):
            mgf.write(spectra_chunk, output_path)
            print(f"{len(spectra_chunk)} spectra with charge written to {output_path}.")
        else:
            print(f"{output_path} already exists.")

        idxs_0 = {idx: spectrum["params"]["title"] for idx, spectrum in enumerate(spectra_chunk)}
        idxs_0 = pd.Series(idxs_0).reset_index()
        idxs_0.columns = ["idx_0", "title"]
        idxs_0["filename"] = output_fname
        spectra_idxs_0.append(idxs_0)

    # Remove the original MGF file after splitting
    # os.remove(input_path)
    # print(f"Original file {input_path} removed.")

spectra_idxs_0 = pd.concat(spectra_idxs_0, axis=0).reset_index(drop=True)

results_df = pd.merge(results_df, spectra_idxs_0, on="title")
results_df["spectrum_id"] = results_df["filename"] + ":" + results_df["idx_0"].astype(str)

results_df["peptide"] = results_df["peptide"].apply(format_peptide_notation)
sequences_true = results_df[["peptide", "spectrum_id"]]
sequences_true = sequences_true.rename({"peptide": "seq"}, axis=1)

labels_path = os.path.join(DATASET_STORAGE_DIR, dset_name, "labels.csv")
sequences_true.to_csv(labels_path, index=False)

# Add dataset tags to the dataset_tags file
collect_dataset_tags(config)

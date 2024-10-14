DATASET_DIR := ./sample_data

# Find all dataset directories
DATASETS := $(shell find $(DATASET_DIR) -mindepth 1 -maxdepth 1 -type d)

# remove all container.sif files
clean:
	find . -type f -name "container.sif" -exec rm -f {} \;

build_apptainer:
	yes | ./build_apptainer_images.sh

# Optionally specify the algorithm to run
# For example:
# make run_all algorithm=instanovo
run_all:
	@for dataset in $(DATASETS); do \
		./run.sh $$dataset $(algorithm); \
	done


streamlit:
	streamlit run dashboard.py --server.enableCORS false --server.enableXsrfProtection false
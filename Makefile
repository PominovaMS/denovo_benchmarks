DATASET_DIR := ./sample_data
RUN_SCRIPT := ./run.sh

# Find all dataset directories
DATASETS := $(shell find $(DATASET_DIR) -mindepth 1 -maxdepth 1 -type d)

build_apptainer:
	yes | ./build_apptainer_images.sh

run:
	$(RUN_SCRIPT) ./sample_data/9_species_human 

instanovo:
	$(RUN_SCRIPT) ./sample_data/9_species_human instanovo

# Define a target to run ./run.sh for each dataset, optionally with an algorithm
# For example:
# make run_all algorithm=instanovo
run_all:
	@for dataset in $(DATASETS); do \
		echo "Running benchmark for dataset: $$dataset"; \
		$(RUN_SCRIPT) $$dataset  $(algorithm); \
	done

streamlit:
	streamlit run dashboard.py --server.enableCORS false --server.enableXsrfProtection false
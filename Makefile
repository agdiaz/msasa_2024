run-all-experiments:
	nextflow run main.nf --input_dir ./datasets

run-experiment:
	# 1. Coincidences
	echo "Coincidences Strategy" > $(INPUT_FASTA_FILE).msasa.timing.log
	date +"%Y-%m-%d %H:%M:%S" >> $(INPUT_FASTA_FILE).msasa.timing.log

	/usr/bin/time -p python msasa_cli.py $(INPUT_FASTA_FILE) $(INPUT_FASTA_FILE).msasa.1.coincidences.aln \
		--log_file $(INPUT_FASTA_FILE).msasa.1.coincidences.log \
		--extend \
		--temp 10 \
		--cooling_rate 0.99 \
		--min_temp 0.001 \
		--max_no_changes 500 \
		--quality_function coincidences 2>> $(INPUT_FASTA_FILE).msasa.timing.log

	# 2. Identity
	echo "Identity Strategy" >> $(INPUT_FASTA_FILE).msasa.timing.log
	date +"%Y-%m-%d %H:%M:%S" >> $(INPUT_FASTA_FILE).msasa.timing.log

	/usr/bin/time -p python msasa_cli.py $(INPUT_FASTA_FILE) $(INPUT_FASTA_FILE).msasa.2.identity.aln \
		--log_file $(INPUT_FASTA_FILE).msasa.2.identity.log \
		--extend \
		--temp 5 \
		--cooling_rate 0.99 \
		--min_temp 0.0001 \
		--max_no_changes 500 \
		--quality_function identity 2>> $(INPUT_FASTA_FILE).msasa.timing.log

	# 3. Similarity
	echo "Similarity Strategy (Blosum62)" >> $(INPUT_FASTA_FILE).msasa.timing.log
	date +"%Y-%m-%d %H:%M:%S" >> $(INPUT_FASTA_FILE).msasa.timing.log

	/usr/bin/time -p python msasa_cli.py $(INPUT_FASTA_FILE) $(INPUT_FASTA_FILE).msasa.3.similarity.aln \
		--log_file $(INPUT_FASTA_FILE).msasa.3.similarity.log \
		--extend \
		--temp 5 \
		--cooling_rate 0.99 \
		--min_temp 0.01 \
		--max_no_changes 500 \
		--quality_function similarity 2>> $(INPUT_FASTA_FILE).msasa.timing.log

	# 4. Global strategy
	echo "Global Strategy" >> $(INPUT_FASTA_FILE).msasa.timing.log
	date +"%Y-%m-%d %H:%M:%S" >> $(INPUT_FASTA_FILE).msasa.timing.log

	/usr/bin/time -p python msasa_cli.py $(INPUT_FASTA_FILE) $(INPUT_FASTA_FILE).msasa.4.global.aln \
		--log_file $(INPUT_FASTA_FILE).msasa.4.global.log \
		--extend \
		--temp 10 \
		--cooling_rate 0.99 \
		--min_temp 0.01 \
		--max_no_changes 500 \
		--quality_function global 2>> $(INPUT_FASTA_FILE).msasa.timing.log

	# 5. Local Strategy
	echo "Local Strategy" >> $(INPUT_FASTA_FILE).msasa.timing.log
	date +"%Y-%m-%d %H:%M:%S" >> $(INPUT_FASTA_FILE).msasa.timing.log

	/usr/bin/time -p python msasa_cli.py $(INPUT_FASTA_FILE) $(INPUT_FASTA_FILE).msasa.5.local.aln \
		--log_file $(INPUT_FASTA_FILE).msasa.5.local.log \
		--extend \
		--temp 10 \
		--cooling_rate 0.99 \
		--min_temp 0.01 \
		--max_no_changes 500 \
		--quality_function local 2>> $(INPUT_FASTA_FILE).msasa.timing.log

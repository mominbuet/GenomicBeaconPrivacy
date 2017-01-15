#Number of queries from adversary
snplist=500 1000 2000 4000 6000 8000 10000 10000 20000 40000 60000 75000 100000 125000 100000 175000 200000 250000 300000
#False positives considered
fplist=0.005 0.006 0.007 0.008 0.009 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.1 0.5 0.75
#deltas
deltalist=1e-06  0.001 0.01
#relatives
philist=0.5 0.25 0.125
#bias for privacy
bias=0.50

path=/import/helium-share/csgrad/azizmma/bustamante/analysis/2016_12_23
#delta 0.000001 0.001 0.01 0.1
#snplist=500 1000 2000 3000 4000 5000 6000 7000 8000 10000 12000 14000 16000 20000 30000 40000
#all: beacon.1000 makeQueries simplepower.png
all: all
#beacon.1000 makeRelatives makeQueryInds makeQueryRelatives makeAltAlleleFiles makeQueries simplepower.png
# beacon.1000 makeRelatives makeQueryInds makeQueryRelatives makeAltAlleleFiles makeQueries simplepower.png

rawdata/Frequencies.txt: ../../scripts/simulate_beacon_saveall_param.r
	mkdir -p $(@D)
	Rscript $< $(@D) 4000000 20000

beacon.2000: ../../scripts/generate_beacon_batches.r $(path)/rawdata/Frequencies.txt
	Rscript $< $@ $(path)/rawdata/Batch[1-9].txt $(path)/rawdata/Batch10.txt $(path)/rawdata/Batch1[1-9].txt $(path)/rawdata/Batch20.txt

beacon.1500: ../../scripts/generate_beacon_batches.r $(path)/rawdata/Frequencies.txt
	Rscript $< $@ $(path)/rawdata/Batch[1-9].txt $(path)/rawdata/Batch10.txt $(path)/rawdata/Batch1[1-5].txt

beacon.1000: ../../scripts/generate_beacon_batches.r $(path)/rawdata/Frequencies.txt
	Rscript $< $@ $(path)/rawdata/Batch[1-9].txt $(path)/rawdata/Batch10.txt

beacon.1200: ../../scripts/generate_beacon_batches.r $(path)/rawdata/Frequencies.txt
	Rscript $< $@ $(path)/rawdata/Batch[1-9].txt $(path)/rawdata/Batch10.txt $(path)/rawdata/Batch1[1-2].txt
	
makeRelatives: ../../scripts/generate_relatives.r $(path)/rawdata/Frequencies.txt
	#ifeq ($(relatives),1)
		for phi in $(philist); do \
			for bno in 1 2 3 4 5 6 7 8 11 12 13 14 15 16 17 18 19 20; do \
				Rscript $< $(word 2, $^) $$phi $(path)/rawdata/Batch$$bno.txt; \
			done ; \
		done ; \
	#endif
	@touch $@

makeQueryInds: ../../scripts/generate_queryinds.r $(path)/rawdata/Frequencies.txt
	for delta in $(deltalist); do \
		for bno in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do \
			Rscript $< $(word 2, $^) $$delta $(path)/rawdata/Batch$$bno.txt; \
	done ; \
	done
	@touch $@

makeQueryRelatives: ../../scripts/generate_queryinds.r $(path)/rawdata/Frequencies.txt
	#ifeq ($(relatives),1)
		for delta in $(deltalist); do \
			for file in rawdata/*related*; do \
				echo $$file $$delta; \
				Rscript $< $(word 2, $^) $$delta $$file; \
			done ; \
		done;\
	#endif
	@touch $@

makeAltAlleleFiles: ../../scripts/create_queryfiles_batches.r 
	for delta in $(deltalist); do \
		mkdir -p delta_$${delta} ; \
		Rscript $< delta_$${delta}/In 1 $(path)/rawdata/Query_$${delta}_Batch1.txt;  \
		Rscript $< delta_$${delta}/In 101 $(path)/rawdata/Query_$${delta}_Batch2.txt; \
		Rscript $< delta_$${delta}/In 201 $(path)/rawdata/Query_$${delta}_Batch3.txt; \
		Rscript $< delta_$${delta}/In 301 $(path)/rawdata/Query_$${delta}_Batch4.txt; \
		Rscript $< delta_$${delta}/In 401 $(path)/rawdata/Query_$${delta}_Batch5.txt; \
		Rscript $< delta_$${delta}/In 501 $(path)/rawdata/Query_$${delta}_Batch6.txt; \
		Rscript $< delta_$${delta}/In 601 $(path)/rawdata/Query_$${delta}_Batch7.txt; \
		Rscript $< delta_$${delta}/In 701 $(path)/rawdata/Query_$${delta}_Batch8.txt; \
		Rscript $< delta_$${delta}/Out 1 $(path)/rawdata/Query_$${delta}_Batch11.txt; \
		Rscript $< delta_$${delta}/Out 101 $(path)/rawdata/Query_$${delta}_Batch12.txt; \
		Rscript $< delta_$${delta}/Out 201 $(path)/rawdata/Query_$${delta}_Batch13.txt; \
		Rscript $< delta_$${delta}/Out 301 $(path)/rawdata/Query_$${delta}_Batch14.txt; \
		Rscript $< delta_$${delta}/Out 401 $(path)/rawdata/Query_$${delta}_Batch15.txt; \
		Rscript $< delta_$${delta}/Out 501 $(path)/rawdata/Query_$${delta}_Batch16.txt; \
		Rscript $< delta_$${delta}/Out 601 $(path)/rawdata/Query_$${delta}_Batch17.txt; \
		Rscript $< delta_$${delta}/Out 701 $(path)/rawdata/Query_$${delta}_Batch18.txt; \
		#
		#	for phi in $(philist); do \
		#		mkdir -p delta_$${delta}/phi_$${phi} ; \
		#		Rscript $< delta_$${delta}/phi_$${phi}/RelIn 1 $(path)/rawdata/Query_$${delta}_Batch1.txt.related_$${phi};  \
		#		Rscript $< delta_$${delta}/phi_$${phi}/RelIn 101 $(path)/rawdata/Query_$${delta}_Batch2.txt.related_$${phi}; \
		#	done; \
		
	done
	@touch $@

makeQueries1000: ../../scripts/findhits.r beacon.1000 makeAltAlleleFiles
	for delta in $(deltalist); do \
		Rscript $< $(word 2,$^) delta_$${delta}/*.altalleles  ; \
		#
		#	for phi in $(philist); do \
		#		Rscript $< $(word 2,$^) delta_$${delta}/phi_$${phi}/*.altalleles  ; \
		#	done; \
		
	done
	@touch $@

getPower1000: ../../scripts/power_altalleles.r makeQueries1000
	
	for nsnps in $(snplist); do \
		echo $$nsnps ; \
		for delta in $(deltalist); do \
			echo "delta" $$delta ; \
			Rscript $< delta_$${delta}/In 400 delta_$${delta}/Out 400 $$nsnps $(fplist) > delta_$${delta}/out.txt.$$nsnps; \
			#
			#	for phi in $(philist); do \
			#		echo "phi" $$phi ; \
			#		Rscript $< delta_$${delta}/phi_$${phi}/RelIn 200 delta_$${delta}/Out 200 $$nsnps $(fplist) ; \
			#	done ; \
			
		done ; \
	done
	@touch $@

allresults1000.dump: getPower1000
		rm -f $@ ;
		echo "Nsnps FPR Power Delta Phi" > $@ ;
		for delta in $(deltalist); do \
			cat delta_$${delta}/out.txt.* | awk -v d=$$delta -v phi=1 '{print $$0,d,phi}' >> $@ ; \
			#
			#	for phi in $(philist); do \
			#		cat delta_$${delta}/phi_$${phi}/out.txt.* | awk -v d=$$delta -v phi=$$phi '{print $$0,d,phi}' >> $@ ; \
			#	done ; \
			
		done ; 
		@touch $@

		
makeRandResponse: ../../scripts/rand_response.r beacon.1000  makeAltAlleleFiles
	for delta in $(deltalist); do \
		Rscript $< $(word 2,$^) $(bias) delta_$${delta}/*.altalleles  ; \
		#
		#	for phi in $(philist); do \
		#		Rscript $< $(word 2,$^) $(bias) delta_$${delta}/phi_$${phi}/*.altalleles  ; \
		#	done; \
		
	done
	@touch $@
	
getPowerRand: ../../scripts/power_altalleles_rand.r makeRandResponse
	
	for nsnps in $(snplist); do \
			echo $$nsnps ; \
			for delta in $(deltalist); do \
				echo "delta" $$delta ; \
				Rscript $< delta_$${delta}/In 400 delta_$${delta}/Out 400 $$nsnps $(fplist) > delta_$${delta}/out_rand.txt.$$nsnps; \
				#
				#	for phi in $(philist); do \
				#		echo "phi" $$phi ; \
				#		Rscript $< delta_$${delta}/phi_$${phi}/RelIn 200 delta_$${delta}/Out 200 $$nsnps $(fplist) > delta_$${delta}/phi_$${phi}/out_rand.txt.$$nsnps; \
				#	done ; \
				
			done ; \
	done
	@touch $@
	
allresultsRand50.dump: getPowerRand
		rm -f $@ ;
		echo "Nsnps FPR Power Delta Phi" > $@ ;
		for delta in $(deltalist); do \
			cat delta_$${delta}/out_rand.txt.* | awk -v d=$$delta -v phi=1 '{print $$0,d,phi}' >> $@ ; \
			#
			#	for phi in $(philist); do \
			#		cat delta_$${delta}/phi_$${phi}/out_rand.txt.* | awk -v d=$$delta -v phi=$$phi '{print $$0,d,phi}' >> $@ ; \
			#	done ; \
			
		done ; 
		@touch $@
		
simplepower.png fp_vs_power.png phi_effect.png delta_effect.png: ../../scripts/newsimulations_plot.r allresults.dump 
	Rscript $^

clean:
	find . -name out.txt.* -exec rm -rf {} \;
	find . -name out_rand.txt.* -exec rm -rf {} \;
	#rm -rf allresults*
	rm -rf getPowerRand
	rm -rf makeRandResponse
	rm -rf *.png


for a in $(seq 1 1 10); do
	for b in $(seq 1 1 10); do
		./GaAs_ED_prism_kinematic.py $a $b
		mv image.png images_for_calibration_wide/image_$a$b.png
		done
	done

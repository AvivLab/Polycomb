cd /polycomb
num=2020_167
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.4 0.5 6.0 0.15 2020_167
python3 upload.py polycomb $num results/$num

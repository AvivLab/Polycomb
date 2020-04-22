cd /polycomb
num=2020_206
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.5 0.5 4.0 0.1 2020_206
python3 upload.py polycomb $num results/$num

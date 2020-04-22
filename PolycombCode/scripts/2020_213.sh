cd /polycomb
num=2020_213
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.4 0.5 1.0 0.1 2020_213
python3 upload.py polycomb $num results/$num

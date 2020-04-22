cd /polycomb
num=2020_192
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.3 0.5 6.0 0.1 2020_192
python3 upload.py polycomb $num results/$num

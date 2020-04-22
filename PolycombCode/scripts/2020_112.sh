cd /polycomb
num=2020_112
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.3 2.0 6.0 0.1 2020_112
python3 upload.py polycomb $num results/$num

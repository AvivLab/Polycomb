cd /polycomb
num=2020_120
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.3 2.0 4.0 0.1 2020_120
python3 upload.py polycomb $num results/$num

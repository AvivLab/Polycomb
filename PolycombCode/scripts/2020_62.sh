cd /polycomb
num=2020_62
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.5 1.0 6.0 0.2 2020_62
python3 upload.py polycomb $num results/$num
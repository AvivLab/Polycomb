cd /polycomb
num=2020_35
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.5 1.0 6.0 0.1 2020_35
python3 upload.py polycomb $num results/$num

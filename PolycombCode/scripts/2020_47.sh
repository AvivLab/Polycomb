cd /polycomb
num=2020_47
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.3 1.0 1.0 0.1 2020_47
python3 upload.py polycomb $num results/$num

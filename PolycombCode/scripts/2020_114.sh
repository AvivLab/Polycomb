cd /polycomb
num=2020_114
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.4 2.0 6.0 0.1 2020_114
python3 upload.py polycomb $num results/$num

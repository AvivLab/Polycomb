cd /polycomb
num=2020_116
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.5 2.0 6.0 0.1 2020_116
python3 upload.py polycomb $num results/$num

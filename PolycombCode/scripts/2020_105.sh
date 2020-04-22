cd /polycomb
num=2020_105
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.4 2.0 1.0 0.15 2020_105
python3 upload.py polycomb $num results/$num

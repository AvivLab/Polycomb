cd /polycomb
num=2020_135
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.5 2.0 1.0 0.1 2020_135
python3 upload.py polycomb $num results/$num

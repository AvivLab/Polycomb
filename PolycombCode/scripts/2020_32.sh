cd /polycomb
num=2020_32
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.4 1.0 6.0 0.1 2020_32
python3 upload.py polycomb $num results/$num

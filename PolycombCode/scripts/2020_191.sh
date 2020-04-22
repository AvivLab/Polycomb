cd /polycomb
num=2020_191
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.3 0.5 6.0 0.1 2020_191
python3 upload.py polycomb $num results/$num

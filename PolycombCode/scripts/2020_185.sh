cd /polycomb
num=2020_185
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.4 0.5 1.0 0.15 2020_185
python3 upload.py polycomb $num results/$num

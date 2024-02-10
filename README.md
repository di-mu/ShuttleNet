# About result data
* Filename format: `res{testing-data-car-id}-i{initial-data-car-id}`
* First column: Reliability (PRR)
* Second column: Normalized Throughput
* Row index: Loop number of shuttle operation

# About code
* To compile
  `gcc shuttle-net.c -o shuttle-net`

* Example usage
  `./shuttle-net initial-data/x380.csv test-data/test380.csv result380.csv`

* Configurable modes and parameters
  `KNN`, `ADR`, `PROB`, `STAT`, `NOGPS`.
  Required Reliability (e.g., `0.8`, `0.9`)

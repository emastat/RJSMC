# Official Simulation Results

This document summarizes the results of running the SMC algorithm with different Update Interval lengths.

## Simulation Setup

- **Number of particles**: 5000
- **Number of iterations**: 50000
- **Burn-in**: 20000
- **Thinning**: 5
- **Method**: turcotte
- **Total observations**: 4683
- **Number of true breakpoints**: 52
- **Number of V states**: 5

## Results Summary

### Performance Summary Table

| length_UI | Run Time (s) | Breakpoint Accuracy | V State Accuracy | Avg Error (Breakpoints) | Avg Error (V State) |
|-----------|---------------|---------------------|------------------|-------------------------|---------------------|
| 0.5 | 884.89 | 0.9895 | 0.9932 | 0.0105 | 0.0068 |
| 1.0 | 465.29 | 0.9922 | 0.9910 | 0.0078 | 0.0090 |
| 2.0 | 196.82 | 0.9960 | 0.9955 | 0.0040 | 0.0045 |
| 5.0 | 112.62 | 0.9968 | 0.9934 | 0.0032 | 0.0066 |

## Run 1: length_UI = 0.5

### Run Configuration

- **Update Interval length**: 0.5
- **Number of Update Intervals**: 95
- **Number of particles**: 5000
- **Number of evaluation intervals**: 4778
- **Number of observations**: 4683
- **Computation time**: 884.89 seconds

### Breakpoint Detection Results

**Confusion Matrix** (rows = true, columns = estimated):

```
       Est_0 Est_1 Est_2
True_0  4691    36     0
True_1    13    37     1
True_2     0     0     0
```

- **Total intervals evaluated**: 4778
- **Correct predictions**: 4728
- **Accuracy**: 0.9895 (98.95%)
- **Average error**: 0.0105

**Interpretation**: The confusion matrix shows how many intervals were correctly/incorrectly classified in terms of breakpoint count (0, 1, or 2 breakpoints per interval).

### V State Classification Results

**Confusion Matrix** (rows = true, columns = estimated):

```
         Est_V_1 Est_V_2 Est_V_3 Est_V_4 Est_V_5
True_V_1     672       6       0       3       0
True_V_2       1     138       0       3       0
True_V_3       0       2      48       0       0
True_V_4       5       9       1    1641       0
True_V_5       0       2       0       0    2152
```

- **Total observations evaluated**: 4683
- **Correct predictions**: 4651
- **Accuracy**: 0.9932 (99.32%)
- **Average error**: 0.0068

**Interpretation**: The confusion matrix shows how many observations were correctly/incorrectly classified in terms of V state (1 to 5).

---

## Run 2: length_UI = 1.0

### Run Configuration

- **Update Interval length**: 1
- **Number of Update Intervals**: 48
- **Number of particles**: 5000
- **Number of evaluation intervals**: 4731
- **Number of observations**: 4683
- **Computation time**: 465.29 seconds

### Breakpoint Detection Results

**Confusion Matrix** (rows = true, columns = estimated):

```
       Est_0 Est_1 Est_2
True_0  4657    26     0
True_1     9    35     1
True_2     0     1     2
```

- **Total intervals evaluated**: 4731
- **Correct predictions**: 4694
- **Accuracy**: 0.9922 (99.22%)
- **Average error**: 0.0078

**Interpretation**: The confusion matrix shows how many intervals were correctly/incorrectly classified in terms of breakpoint count (0, 1, or 2 breakpoints per interval).

### V State Classification Results

**Confusion Matrix** (rows = true, columns = estimated):

```
         Est_V_1 Est_V_2 Est_V_3 Est_V_4 Est_V_5
True_V_1     675       6       0       0       0
True_V_2       0     135       4       3       0
True_V_3       0       0      44       6       0
True_V_4       2       9       3    1642       0
True_V_5       0       5       0       4    2145
```

- **Total observations evaluated**: 4683
- **Correct predictions**: 4641
- **Accuracy**: 0.9910 (99.10%)
- **Average error**: 0.0090

**Interpretation**: The confusion matrix shows how many observations were correctly/incorrectly classified in terms of V state (1 to 5).

---

## Run 3: length_UI = 2.0

### Run Configuration

- **Update Interval length**: 2
- **Number of Update Intervals**: 24
- **Number of particles**: 5000
- **Number of evaluation intervals**: 4707
- **Number of observations**: 4683
- **Computation time**: 196.82 seconds

### Breakpoint Detection Results

**Confusion Matrix** (rows = true, columns = estimated):

```
       Est_0 Est_1 Est_2
True_0  4651    10     0
True_1     7    33     1
True_2     0     1     4
```

- **Total intervals evaluated**: 4707
- **Correct predictions**: 4688
- **Accuracy**: 0.9960 (99.60%)
- **Average error**: 0.0040

**Interpretation**: The confusion matrix shows how many intervals were correctly/incorrectly classified in terms of breakpoint count (0, 1, or 2 breakpoints per interval).

### V State Classification Results

**Confusion Matrix** (rows = true, columns = estimated):

```
         Est_V_1 Est_V_2 Est_V_3 Est_V_4 Est_V_5
True_V_1     680       1       0       0       0
True_V_2       0     139       0       3       0
True_V_3       0       0      50       0       0
True_V_4       9       3       0    1644       0
True_V_5       0       2       0       3    2149
```

- **Total observations evaluated**: 4683
- **Correct predictions**: 4662
- **Accuracy**: 0.9955 (99.55%)
- **Average error**: 0.0045

**Interpretation**: The confusion matrix shows how many observations were correctly/incorrectly classified in terms of V state (1 to 5).

---

## Run 4: length_UI = 5.0

### Run Configuration

- **Update Interval length**: 5
- **Number of Update Intervals**: 10
- **Number of particles**: 5000
- **Number of evaluation intervals**: 4693
- **Number of observations**: 4683
- **Computation time**: 112.62 seconds

### Breakpoint Detection Results

**Confusion Matrix** (rows = true, columns = estimated):

```
       Est_0 Est_1 Est_2
True_0  4645     7     0
True_1     8    23     0
True_2     0     0    10
```

- **Total intervals evaluated**: 4693
- **Correct predictions**: 4678
- **Accuracy**: 0.9968 (99.68%)
- **Average error**: 0.0032

**Interpretation**: The confusion matrix shows how many intervals were correctly/incorrectly classified in terms of breakpoint count (0, 1, or 2 breakpoints per interval).

### V State Classification Results

**Confusion Matrix** (rows = true, columns = estimated):

```
         Est_V_1 Est_V_2 Est_V_3 Est_V_4 Est_V_5
True_V_1     668       4       0       9       0
True_V_2       0     139       0       3       0
True_V_3       0       0      44       6       0
True_V_4       1       6       0    1649       0
True_V_5       0       2       0       0    2152
```

- **Total observations evaluated**: 4683
- **Correct predictions**: 4652
- **Accuracy**: 0.9934 (99.34%)
- **Average error**: 0.0066

**Interpretation**: The confusion matrix shows how many observations were correctly/incorrectly classified in terms of V state (1 to 5).

---

## Comparative Analysis

### Breakpoint Detection Performance

| length_UI | Accuracy | Avg Error |
|-----------|----------|-----------|
| 0.5 | 0.9895 | 0.0105 |
| 1.0 | 0.9922 | 0.0078 |
| 2.0 | 0.9960 | 0.0040 |
| 5.0 | 0.9968 | 0.0032 |

### V State Classification Performance

| length_UI | Accuracy | Avg Error |
|-----------|----------|-----------|
| 0.5 | 0.9932 | 0.0068 |
| 1.0 | 0.9910 | 0.0090 |
| 2.0 | 0.9955 | 0.0045 |
| 5.0 | 0.9934 | 0.0066 |

### Computational Performance

| length_UI | Run Time (seconds) | Update Intervals |
|-----------|---------------------|------------------|
| 0.5 | 884.89 | 95 |
| 1.0 | 465.29 | 48 |
| 2.0 | 196.82 | 24 |
| 5.0 | 112.62 | 10 |

## Notes

- **Breakpoint Confusion Matrix**: 3x3 matrix (rows = true, columns = estimated breakpoint count 0/1/2).
- **V State Confusion Matrix**: 5x5 matrix (rows = true, columns = estimated V state).
- **Average Error**: Proportion of incorrect predictions.



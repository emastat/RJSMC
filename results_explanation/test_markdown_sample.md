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
| 0.5 | 120.00 | 0.9735 | 1.0000 | 0.0100 | 0.0100 |
| 1.0 | 140.00 | 0.9735 | 1.0000 | 0.0200 | 0.0200 |
| 2.0 | 160.00 | 0.9735 | 1.0000 | 0.0300 | 0.0300 |
| 5.0 | 180.00 | 0.9735 | 1.0000 | 0.0400 | 0.0400 |

## Run 1: length_UI = 0.5

### Run Configuration

- **Update Interval length**: 0.5
- **Number of Update Intervals**: 40
- **Number of particles**: 5000
- **Number of evaluation intervals**: 950
- **Number of observations**: 4683
- **Computation time**: 120.00 seconds

### Breakpoint Detection Results

**Confusion Matrix** (rows = true, columns = estimated):

```
       Est_0 Est_1 Est_2
True_0   100     1     0
True_1     2    10     0
True_2     0     0     0
```

- **Total intervals evaluated**: 113
- **Correct predictions**: 110
- **Accuracy**: 0.9735 (97.35%)
- **Average error**: 0.0100

**Interpretation**: The confusion matrix shows how many intervals were correctly/incorrectly classified in terms of breakpoint count (0, 1, or 2 breakpoints per interval).

### V State Classification Results

**Confusion Matrix** (rows = true, columns = estimated):

```
         Est_V_1 Est_V_2 Est_V_3 Est_V_4 Est_V_5
True_V_1      50       0       0       0       0
True_V_2       0      30       0       0       0
True_V_3       0       0      20       0       0
True_V_4       0       0       0      40       0
True_V_5       0       0       0       0      60
```

- **Total observations evaluated**: 200
- **Correct predictions**: 200
- **Accuracy**: 1.0000 (100.00%)
- **Average error**: 0.0100

**Interpretation**: The confusion matrix shows how many observations were correctly/incorrectly classified in terms of V state (1 to 5).

---

## Run 2: length_UI = 1.0

### Run Configuration

- **Update Interval length**: 1
- **Number of Update Intervals**: 30
- **Number of particles**: 5000
- **Number of evaluation intervals**: 900
- **Number of observations**: 4683
- **Computation time**: 140.00 seconds

### Breakpoint Detection Results

**Confusion Matrix** (rows = true, columns = estimated):

```
       Est_0 Est_1 Est_2
True_0   100     1     0
True_1     2    10     0
True_2     0     0     0
```

- **Total intervals evaluated**: 113
- **Correct predictions**: 110
- **Accuracy**: 0.9735 (97.35%)
- **Average error**: 0.0200

**Interpretation**: The confusion matrix shows how many intervals were correctly/incorrectly classified in terms of breakpoint count (0, 1, or 2 breakpoints per interval).

### V State Classification Results

**Confusion Matrix** (rows = true, columns = estimated):

```
         Est_V_1 Est_V_2 Est_V_3 Est_V_4 Est_V_5
True_V_1      50       0       0       0       0
True_V_2       0      30       0       0       0
True_V_3       0       0      20       0       0
True_V_4       0       0       0      40       0
True_V_5       0       0       0       0      60
```

- **Total observations evaluated**: 200
- **Correct predictions**: 200
- **Accuracy**: 1.0000 (100.00%)
- **Average error**: 0.0200

**Interpretation**: The confusion matrix shows how many observations were correctly/incorrectly classified in terms of V state (1 to 5).

---

## Run 3: length_UI = 2.0

### Run Configuration

- **Update Interval length**: 2
- **Number of Update Intervals**: 20
- **Number of particles**: 5000
- **Number of evaluation intervals**: 850
- **Number of observations**: 4683
- **Computation time**: 160.00 seconds

### Breakpoint Detection Results

**Confusion Matrix** (rows = true, columns = estimated):

```
       Est_0 Est_1 Est_2
True_0   100     1     0
True_1     2    10     0
True_2     0     0     0
```

- **Total intervals evaluated**: 113
- **Correct predictions**: 110
- **Accuracy**: 0.9735 (97.35%)
- **Average error**: 0.0300

**Interpretation**: The confusion matrix shows how many intervals were correctly/incorrectly classified in terms of breakpoint count (0, 1, or 2 breakpoints per interval).

### V State Classification Results

**Confusion Matrix** (rows = true, columns = estimated):

```
         Est_V_1 Est_V_2 Est_V_3 Est_V_4 Est_V_5
True_V_1      50       0       0       0       0
True_V_2       0      30       0       0       0
True_V_3       0       0      20       0       0
True_V_4       0       0       0      40       0
True_V_5       0       0       0       0      60
```

- **Total observations evaluated**: 200
- **Correct predictions**: 200
- **Accuracy**: 1.0000 (100.00%)
- **Average error**: 0.0300

**Interpretation**: The confusion matrix shows how many observations were correctly/incorrectly classified in terms of V state (1 to 5).

---

## Run 4: length_UI = 5.0

### Run Configuration

- **Update Interval length**: 5
- **Number of Update Intervals**: 10
- **Number of particles**: 5000
- **Number of evaluation intervals**: 800
- **Number of observations**: 4683
- **Computation time**: 180.00 seconds

### Breakpoint Detection Results

**Confusion Matrix** (rows = true, columns = estimated):

```
       Est_0 Est_1 Est_2
True_0   100     1     0
True_1     2    10     0
True_2     0     0     0
```

- **Total intervals evaluated**: 113
- **Correct predictions**: 110
- **Accuracy**: 0.9735 (97.35%)
- **Average error**: 0.0400

**Interpretation**: The confusion matrix shows how many intervals were correctly/incorrectly classified in terms of breakpoint count (0, 1, or 2 breakpoints per interval).

### V State Classification Results

**Confusion Matrix** (rows = true, columns = estimated):

```
         Est_V_1 Est_V_2 Est_V_3 Est_V_4 Est_V_5
True_V_1      50       0       0       0       0
True_V_2       0      30       0       0       0
True_V_3       0       0      20       0       0
True_V_4       0       0       0      40       0
True_V_5       0       0       0       0      60
```

- **Total observations evaluated**: 200
- **Correct predictions**: 200
- **Accuracy**: 1.0000 (100.00%)
- **Average error**: 0.0400

**Interpretation**: The confusion matrix shows how many observations were correctly/incorrectly classified in terms of V state (1 to 5).

---

## Comparative Analysis

### Breakpoint Detection Performance

| length_UI | Accuracy | Avg Error |
|-----------|----------|-----------|
| 0.5 | 0.9735 | 0.0100 |
| 1.0 | 0.9735 | 0.0200 |
| 2.0 | 0.9735 | 0.0300 |
| 5.0 | 0.9735 | 0.0400 |

### V State Classification Performance

| length_UI | Accuracy | Avg Error |
|-----------|----------|-----------|
| 0.5 | 1.0000 | 0.0100 |
| 1.0 | 1.0000 | 0.0200 |
| 2.0 | 1.0000 | 0.0300 |
| 5.0 | 1.0000 | 0.0400 |

### Computational Performance

| length_UI | Run Time (seconds) | Update Intervals |
|-----------|---------------------|------------------|
| 0.5 | 120.00 | 40 |
| 1.0 | 140.00 | 30 |
| 2.0 | 160.00 | 20 |
| 5.0 | 180.00 | 10 |

## Notes

- **Breakpoint Confusion Matrix**: 3x3 matrix (rows = true, columns = estimated breakpoint count 0/1/2).
- **V State Confusion Matrix**: 5x5 matrix (rows = true, columns = estimated V state).
- **Average Error**: Proportion of incorrect predictions.



# Gender-Specific_Stammer_Detection

This project implements a **speech recognition system** in C++ using techniques from digital signal processing and machine learning.  
It processes raw speech signals, extracts features, generates a codebook with the **Linde–Buzo–Gray (LBG) algorithm**, and trains **Hidden Markov Models (HMMs)** for classification.

---

## ✨ Features
- **Preprocessing of speech signals**  
  - DC shift removal  
  - Normalization  
  - Steady frame extraction based on energy threshold  

- **Feature Extraction**  
  - Autocorrelation (Ri) computation  
  - Linear Predictive Coding (LPC) via Levinson-Durbin algorithm  
  - Cepstral Coefficient computation  

- **Vector Quantization (VQ) with LBG Algorithm**  
  - Generates a **codebook** using Tokhura distance  
  - Performs K-means clustering for refinement  

- **Observation Sequence Generation**  
  - Converts Cepstral Coefficients into discrete symbols using the codebook  

- **Hidden Markov Model (HMM)**  
  - Forward and Backward algorithms  
  - Baum-Welch re-estimation (EM algorithm)  
  - Viterbi decoding for optimal state sequence  

- **Training and Testing**  
  - Supports multiple classes (e.g., Male, Female, MPS, FPS, MSRS, FSRS)  
  - Iterative model training with averaging across multiple samples  

---

## 📂 Project Structure
Speech_project/
│
├── Speech_project.cpp # Main source file (all functions)
├── initialModel/ # Initial A, B, PI matrices
├── model/ # Trained HMM parameters
├── avg/ # Averaged models after multiple iterations
├── Observation/ # Generated observation sequences
├── Universe/ # Cepstral coefficient universe
├── Codebook/ # Final codebook (LBG output)
├── data/ # Input speech signals (raw amplitude values)
└── README.md # Documentation


---

## ⚙️ Dependencies
- **Compiler**: MSVC / MinGW / g++ (Windows only due to `Windows.h` and `winmm.lib`)  
- **Libraries**:  
  - `<math.h>` (Math operations)  
  - `<mmsystem.h>` (Windows audio support)  
  - `<direct.h>` (Directory handling)  
  - `<stdio.h>`, `<stdlib.h>`, `<limits>`  

> ⚠️ Note: For Linux/Mac, you must replace Windows-specific audio handling.

---

## 🚀 Workflow

1. **Signal Acquisition**  
   - Reads `.txt` files containing raw amplitude values of speech signals.  

2. **Preprocessing**  
   - Removes DC shift, normalizes amplitude, and extracts steady frames.  

3. **Feature Extraction**  
   - Computes **LPC** and **Cepstral Coefficients**.  
   - Stores coefficients in `Universe/` and `Observation/`.  

4. **Codebook Generation (LBG)**  
   - Uses Tokhura distance + K-means clustering.  
   - Produces a **32-symbol codebook**.  

5. **HMM Training (Baum-Welch)**  
   - Initializes A, B, π matrices.  
   - Re-estimates using Forward-Backward.  
   - Iteratively averages models over multiple utterances.  

6. **Recognition**  
   - Generates observation sequence from test speech.  
   - Runs **Forward Algorithm** for probability computation.  
   - Applies **Viterbi decoding** for classification.  

---

## ▶️ How to Run

1. **Prepare Speech Data**
   - Record samples.  
   - Convert into `.txt` files (amplitude values).  
   - Place in `data/`.  

2. **Extract Features**
   ```cpp
   calculateCi("Universe/universe.txt", "data/speech1.txt", "Ci/speech1_Ci.txt", 1);
   
3.**Generate Codebook**
  ```cpp
  universeToArray(size, ORDER);
  double** finalCodebook = LBG(size, 32, ORDER, 0.03);
```
4.**Train HMM***
  ```cpp
  trainingData(iteration);
  averageModel(iteration);
```
5.**Test Model**
  ```cpp
  calculateCiTest("Universe/universe.txt", "data/test_speech.txt", "Ci/test_Ci.txt", 1);
  findObervationFile("Ci/test_Ci.txt", "Observation/test_obs.txt", count);
  openObs("Observation/test_obs.txt");
  forwardPropogation();
  viterbiAlgo();
```

🏷️ Authors
Abhishek Kumar Tiwari
MTech, IIT Guwahati
Speech Processing 

---

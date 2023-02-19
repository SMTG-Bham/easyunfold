#


## Wavecar
```python 
Wavecar(
   fnm = 'WAVECAR', lsorbit = False, lgamma = False, gamma_half = 'x',
   omp_num_threads = 1
)
```


---
Class for processing VASP Pseudowavefunction stored in WAVECAR.

This class is a trimmed-down from that of [PyVaspwfc](https://github.com/QijingZheng/VaspBandUnfolding)
by QijingZheng <zqj.kaka@gmail.com>.

Only the functions needed for band unfolding remains.
The original class has many other functionalities such as ELF calculation.

The format of VASP WAVECAR, as shown in
http://www.andrew.cmu.edu/user/feenstra/wavetrans/
---
is:
    Record-length #spin components RTAG(a value specifying the precision)
    #k-points #bands ENCUT(maximum energy for plane waves)
    LatVec-A
    LatVec-B
    LatVec-C
    Loop over spin
       Loop over k-points
          #plane waves, k vector
          Loop over bands
             band energy, band occupation
          End loop over bands
          Loop over bands
             Loop over plane waves
                Plane-wave coefficient
             End loop over plane waves
          End loop over bands
       End loop over k-points
    End loop over spin


**Methods:**


### .isSocWfc
```python
.isSocWfc()
```

---
Is the WAVECAR from an SOC calculation?

### .isGammaWfc
```python
.isGammaWfc()
```

---
Is the WAVECAR from an SOC calculation?

### .readWFHeader
```python
.readWFHeader()
```

---
Read the system information from WAVECAR, which is written in the first
two record.

rec1: recl, nspin, rtag
rec2: nkpts, nbands, encut, ((cell(i,j) i=1, 3), j=1, 3)

### .setWFPrec
```python
.setWFPrec()
```

---
Set wavefunction coefficients precision:
TAG = 45200: single precision complex, np.complex64, or complex(qs)
TAG = 45210: double precision complex, np.complex128, or complex(q)

### .readWFBand
```python
.readWFBand()
```

---
Extract KS energies and Fermi occupations from WAVECAR.

### .gvectors
```python
.gvectors(
   ikpt = 1, force_Gamma = False, check_consistency = True
)
```

---
Generate the G-vectors that satisfies the following relation
(G + k)**2 / 2 < ENCUT

### .readBandCoeff
```python
.readBandCoeff(
   ispin = 1, ikpt = 1, iband = 1, norm = False
)
```

---
Read the planewave coefficients of specified KS states.

### .whereRec
```python
.whereRec(
   ispin = 1, ikpt = 1, iband = 1
)
```

---
Return the rec position for specified KS state.

### .checkIndex
```python
.checkIndex(
   ispin, ikpt, iband
)
```

---
Check if the index is valid!

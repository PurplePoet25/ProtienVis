# 🧬 DNA to Protein Visualizer 🧬

✨ **DNA to Protein Visualizer** is a Python tool that lets you input a DNA sequence, validates it, translates it into amino acids, and visualizes the protein structure in a fun, detailed skeletal view using matplotlib.

## 🌟 Features

* ✅ DNA sequence validation (checks for valid nucleotides, start codon, stop codon)
* ✅ Translation from DNA codons to amino acids
* ✅ Detailed amino acid visualization including side chains (R-groups)
* ✅ Handles most standard amino acids (GLY, ALA, CYS, MET, PHE, LEU, ILE, VAL, SER, THR, ASN, GLN, TYR, TRP, ASP, GLU, HIS, LYS, ARG, PRO)
* 🖥️ GUI input via PyAutoGUI or optional terminal input

## 🛠️ Installation & Setup

1. Clone the repository:

```bash
git clone https://github.com/PurplePoet25/ProtienVis.git
cd ProtienVis
```

2. Install required Python packages:

```bash
pip install matplotlib pyautogui
```

> PyAutoGUI is used for GUI prompts, but you can switch to `input()` for terminal-based input.

## 🎮 How to Use

Run the main script:

```bash
python ProtienVis.py
```

* Enter a DNA sequence (e.g., `ATGGCTTGTGGGTTTAA`).
* The sequence is validated for correct nucleotides, presence of a start codon (`ATG`), and at least one stop codon (`TAA`, `TAG`, `TGA`).
* The DNA is split into codons, translated into amino acids, and visualized as a 2D skeletal diagram with R-groups.

## 🧪 Supported Amino Acids

| Amino Acid    | Abbreviation |
| ------------- | ------------ |
| Glycine       | GLY          |
| Alanine       | ALA          |
| Cysteine      | CYS          |
| Methionine    | MET          |
| Phenylalanine | PHE          |
| Leucine       | LEU          |
| Isoleucine    | ILE          |
| Valine        | VAL          |
| Serine        | SER          |
| Threonine     | THR          |
| Asparagine    | ASN          |
| Glutamine     | GLN          |
| Tyrosine      | TYR          |
| Tryptophan    | TRP          |
| Aspartic Acid | ASP          |
| Glutamic Acid | GLU          |
| Histidine     | HIS          |
| Lysine        | LYS          |
| Arginine      | ARG          |
| Proline       | PRO          |
| Stop Codons   | STOP         |

## 📸 Screenshots

<img width="466" height="195" alt="image" src="https://github.com/user-attachments/assets/d70d0f3e-6551-41d6-bba1-a3eaf411b32c" />
<img width="1242" height="830" alt="image" src="https://github.com/user-attachments/assets/491051f0-bebb-4159-9638-c303058cee6e" />
<img width="1246" height="831" alt="image" src="https://github.com/user-attachments/assets/56166afd-2f05-4a3a-aad2-6bc9f43b67fd" />


## 📁 File Structure

```
ProtienVis/
│
├── ProtienVis.py              # Main Python script
├── README.md            # This documentation

```

## 📝 Example

```python
Enter DNA Sequence: ATGGCTTGTGGGTTTAA
```

This generates a matplotlib figure showing amino acids as colored circles with detailed R-groups.

## ⚠️ Notes

* Long sequences may require larger figure size or horizontal scrolling.
* Unknown codons will appear as `???`.
* PROline is handled specially due to its cyclic structure.

## 🤝 Contributing

Ideas to add more features, improved visualization, or other enhancements are welcome! Open an issue or submit a PR.

## 📜 License

All rights reserved by the author. 🛡️ No copying, modification, or redistribution without permission.

## 👤 Author

Your Name - [GitHub Profile](https://github.com/PurplePoet25)


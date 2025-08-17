import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import pyautogui

CODON_TABLE = {
    "ATA": "ILE", "ATC": "ILE", "ATT": "ILE",
    "ATG": "MET",
    
    "ACA": "THR", "ACC": "THR", "ACG": "THR", "ACT": "THR",
    "AAC": "ASN", "AAT": "ASN",
    "AAA": "LYS", "AAG": "LYS",
    "AGC": "SER", "AGT": "SER",
    "AGA": "ARG", "AGG": "ARG",

    "CTA": "LEU", "CTC": "LEU", "CTG": "LEU", "CTT": "LEU",
    "CCA": "PRO", "CCC": "PRO", "CCG": "PRO", "CCT": "PRO",
    "CAC": "HIS", "CAT": "HIS",
    "CAA": "GLN", "CAG": "GLN",
    "CGA": "ARG", "CGC": "ARG", "CGG": "ARG", "CGT": "ARG",

    "GTA": "VAL", "GTC": "VAL", "GTG": "VAL", "GTT": "VAL",
    "GCA": "ALA", "GCC": "ALA", "GCG": "ALA", "GCT": "ALA",
    "GAC": "ASP", "GAT": "ASP",
    "GAA": "GLU", "GAG": "GLU",
    "GGA": "GLY", "GGC": "GLY", "GGG": "GLY", "GGT": "GLY",

    "TCA": "SER", "TCC": "SER", "TCG": "SER", "TCT": "SER",
    "TTC": "PHE", "TTT": "PHE",
    "TTA": "LEU", "TTG": "LEU",
    "TAC": "TYR", "TAT": "TYR",
    "TAA": "STOP", "TAG": "STOP", "TGA": "STOP",
    "TGC": "CYS", "TGT": "CYS",
    "TGG": "TRP"
}


AA_COLORS = {
    'ALA': '#B0B0B0',
    'ARG': '#00CED1',
    'ASN': '#87CEEB',
    'ASP': '#4169E1',
    'CYS': '#FFD700',
    'GLU': '#00008B',
    'GLN': '#66CDAA',
    'GLY': '#FFFFFF',
    'HIS': '#9370DB',
    'ILE': '#6B8E23',
    'LEU': '#228B22',
    'LYS': '#FF69B4',
    'MET': '#DAA520',
    'PHE': "#B7BC47",
    'PRO': '#D8BFD8',
    'SER': '#B0E0E6',
    'THR': '#00FFFF',
    'TRP': '#4B0082',
    'TYR': '#A52A2A',
    'VAL': '#556B2F'
}


ATOM_COLORS = {
    "C": "black",
    "O": "red",
    "O⁻": "red",
    "N": "blue",
    "N⁺": "blue",
    "H": "gray",
    "S": "yellow"
}

def draw_atom(ax, symbol, x, y):
    face = ATOM_COLORS.get(symbol, 'gray')
    text_color = 'black' if symbol == 'S' else 'white'
    ax.add_patch(patches.Circle((x, y), 0.15, facecolor=face, edgecolor='black'))
    ax.text(x, y, symbol, fontsize=8, ha='center', va='center', color=text_color, weight='bold')

def draw_bond(ax, x1, y1, x2, y2):
    ax.plot([x1, x2], [y1, y2], color='gray', linewidth=1)

def draw_gly(ax, x, y):
    draw_atom(ax, "H", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)

def draw_ala(ax, x, y):
    draw_atom(ax, "C", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)
    draw_atom(ax, "H", x - 0.3, y - 0.9)
    draw_bond(ax, x, y - 0.5, x - 0.3, y - 0.9)
    draw_atom(ax, "H", x + 0.3, y - 0.9)
    draw_bond(ax, x, y - 0.5, x + 0.3, y - 0.9)
    draw_atom(ax, "H", x, y - 1.3)
    draw_bond(ax, x, y - 0.5, x, y - 1.3)

def draw_cys(ax, x, y):
    draw_atom(ax, "C", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)
    draw_atom(ax, "H", x - 0.3, y - 0.9)
    draw_bond(ax, x, y - 0.5, x - 0.3, y - 0.9)
    draw_atom(ax, "H", x + 0.3, y - 0.9)
    draw_bond(ax, x, y - 0.5, x + 0.3, y - 0.9)
    draw_atom(ax, "S", x, y - 1.3)
    draw_bond(ax, x, y - 0.5, x, y - 1.3)
    draw_atom(ax, "H", x + 0.3, y - 1.7)
    draw_bond(ax, x, y - 1.3, x + 0.3, y - 1.7)

def draw_met(ax, x, y):
    atoms = {
        "C1": (0, -0.5), "H1": (-0.3, -0.9), "H2": (0.3, -0.9),
        "C2": (0, -1.3), "H3": (-0.3, -1.7), "H4": (0.3, -1.7),
        "S": (0, -2.1), "C3": (0.4, -2.5),
        "H5": (0.2, -2.9), "H6": (0.6, -2.9), "H7": (0.8, -2.6)
    }
    bonds = [
        ("C1", "H1"), ("C1", "H2"), ("C1", "C2"),
        ("C2", "H3"), ("C2", "H4"), ("C2", "S"),
        ("S", "C3"), ("C3", "H5"), ("C3", "H6"), ("C3", "H7")
    ]
    for label, (dx, dy) in atoms.items():
        draw_atom(ax, label[0], x + dx, y + dy)
    for a1, a2 in bonds:
        x1, y1 = atoms[a1]
        x2, y2 = atoms[a2]
        draw_bond(ax, x + x1, y + y1, x + x2, y + y2)

def draw_phe(ax, x, y):
    # Center Carbon (attached to backbone)
    draw_atom(ax, "C", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)

    # Middle connector Carbon
    draw_atom(ax, "C", x, y - 1.3)
    draw_bond(ax, x, y - 0.5, x, y - 1.3)

    # Hydrogens on central carbon (approx)
    draw_atom(ax, "H", x - 0.4, y - 0.8)
    draw_bond(ax, x, y - 0.5, x - 0.4, y - 0.8)

    draw_atom(ax, "H", x + 0.4, y - 0.8)
    draw_bond(ax, x, y - 0.5, x + 0.4, y - 0.8)

    # Phenyl Ring Positions
    ring = [
        (x - 0.6, y - 1.6),  # left top
        (x - 0.6, y - 2.2),  # left bottom
        (x, y - 2.5),        # bottom
        (x + 0.6, y - 2.2),  # right bottom
        (x + 0.6, y - 1.6),  # right top
        (x, y - 1.3)         # top (shared with middle)
    ]

    # Hydrogens for ring (positions mimicking your primitive data)
    ring_H = [
        (x - 0.75, y - 1.55),  # left top H
        (x - 0.75, y - 2.3),   # left bottom H
        (x, y - 2.7),          # bottom H
        (x + 0.75, y - 2.3),   # right bottom H
        (x + 0.75, y - 1.55)   # right top H
    ]

    # Draw ring bonds & atoms
    for i in range(len(ring)):
        x1, y1 = ring[i]
        x2, y2 = ring[(i + 1) % len(ring)]
        draw_bond(ax, x1, y1, x2, y2)
        if i != 5:  # skip drawing the shared center again
            draw_atom(ax, "C", x1, y1)

    # Draw hydrogens and their bonds
    for (hx, hy), (cx, cy) in zip(ring_H, ring[:-1]):
        draw_atom(ax, "H", hx, hy)
        draw_bond(ax, hx, hy, cx, cy)

def draw_ile(ax, x, y):
    atoms = {
        "C1": (0, -0.5),
        "C2": (0, -1.2),
        "C3": (-0.8, -1.8),   # wide left branch
        "C4": (0, -2.3),      # vertical CH3
        "H1": (0.7, -1.7),    # right-side H on C2

        # Hydrogens on C3 (left CH3)
        "H2": (-1.0, -2.2),
        "H3": (-0.6, -2.2),
        "H4": (-0.8, -2.6),

        # Hydrogens on C4 (bottom CH3)
        "H5": (-0.3, -2.8),
        "H6": (0.3, -2.8),
        "H7": (0, -3.1),
    }

    bonds = [
        ("C1", "C2"), ("C2", "C3"), ("C2", "C4"), ("C2", "H1"),
        ("C3", "H2"), ("C3", "H3"), ("C3", "H4"),
        ("C4", "H5"), ("C4", "H6"), ("C4", "H7")
    ]

    for label, (dx, dy) in atoms.items():
        draw_atom(ax, label[0], x + dx, y + dy)
    for a1, a2 in bonds:
        x1, y1 = atoms[a1]
        x2, y2 = atoms[a2]
        draw_bond(ax, x + x1, y + y1, x + x2, y + y2)


def draw_leu(ax, x, y):
    atoms = {
        "C1": (0, -0.5),
        "C2": (0, -1.2),
        "C3": (-0.9, -1.9),  # wide left CH3
        "C4": (0.9, -1.9),   # wide right CH3
        "H1": (0, -2.4),     # downward H on C2

        # Hydrogens on C3
        "H2": (-1.1, -2.3),
        "H3": (-0.7, -2.3),
        "H4": (-0.9, -2.6),

        # Hydrogens on C4
        "H5": (0.7, -2.3),
        "H6": (1.1, -2.3),
        "H7": (0.9, -2.6),
    }

    bonds = [
        ("C1", "C2"), ("C2", "C3"), ("C2", "C4"), ("C2", "H1"),
        ("C3", "H2"), ("C3", "H3"), ("C3", "H4"),
        ("C4", "H5"), ("C4", "H6"), ("C4", "H7")
    ]

    for label, (dx, dy) in atoms.items():
        draw_atom(ax, label[0], x + dx, y + dy)
    for a1, a2 in bonds:
        x1, y1 = atoms[a1]
        x2, y2 = atoms[a2]
        draw_bond(ax, x + x1, y + y1, x + x2, y + y2)


def draw_val(ax, x, y):
    draw_atom(ax, "C", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)

    # Middle Carbon
    center = (x, y - 1.2)
    draw_atom(ax, "C", *center)
    draw_bond(ax, x, y - 0.5, *center)

    # Side branches
    left = (x - 0.7, y - 1.8)
    right = (x + 0.7, y - 1.8)
    draw_atom(ax, "C", *left)
    draw_atom(ax, "C", *right)
    draw_bond(ax, *center, *left)
    draw_bond(ax, *center, *right)

    # Central H
    draw_atom(ax, "H", x, y - 1.7)
    draw_bond(ax, *center, x, y - 1.7)

    # Fan out Hydrogens - wider and middle one lower
    for i, dx in enumerate([-0.35, 0.0, 0.35]):
        dy = -0.35 if i != 1 else -0.45
        draw_atom(ax, "H", left[0] + dx, left[1] + dy)
        draw_bond(ax, *left, left[0] + dx, left[1] + dy)

        draw_atom(ax, "H", right[0] + dx, right[1] + dy)
        draw_bond(ax, *right, right[0] + dx, right[1] + dy)

def draw_ser(ax, x, y):
    # First Carbon (Cα)
    draw_atom(ax, "C", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)

    # OH Group on side chain
    o = (x, y - 1.2)
    draw_atom(ax, "O", *o)
    draw_bond(ax, x, y - 0.5, *o)

    # Hydroxyl H (from OH)
    draw_atom(ax, "H", x - 0.35, y - 1.5)
    draw_bond(ax, *o, x - 0.35, y - 1.5)

    # First hydrogen on Cα (right side)
    draw_atom(ax, "H", x + 0.5, y - 0.9)
    draw_bond(ax, x, y - 0.5, x + 0.5, y - 0.9)

    # Second hydrogen on Cα (left side)
    draw_atom(ax, "H", x - 0.5, y - 0.9)
    draw_bond(ax, x, y - 0.5, x - 0.5, y - 0.9)

def draw_thr(ax, x, y):
    draw_atom(ax, "C", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)

    # Central Carbon
    center = (x, y - 1.2)
    draw_atom(ax, "C", *center)
    draw_bond(ax, x, y - 0.5, *center)

    # OH Branch
    o = (x, y - 1.9)
    draw_atom(ax, "O", *o)
    draw_bond(ax, *center, *o)
    draw_atom(ax, "H", x - 0.4, y - 2.2)
    draw_bond(ax, *o, x - 0.4, y - 2.2)

    # Methyl branch
    c = (x + 0.8, y - 1.4)
    draw_atom(ax, "C", *c)
    draw_bond(ax, *center, *c)

    # Wider fan of hydrogens with middle lower
    for i, dx in enumerate([-0.35, 0.0, 0.35]):
        dy = -0.3 if i != 1 else -0.45
        draw_atom(ax, "H", c[0] + dx, c[1] + dy)
        draw_bond(ax, *c, c[0] + dx, c[1] + dy)

def draw_asn(ax, x, y):
    # Top carbon (from ASN node)
    draw_atom(ax, "C", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)

    # Two hydrogens on top carbon (corrected)
    draw_atom(ax, "H", x - 0.4, y - 0.95)
    draw_bond(ax, x, y - 0.5, x - 0.4, y - 0.95)

    draw_atom(ax, "H", x + 0.4, y - 0.95)
    draw_bond(ax, x, y - 0.5, x + 0.4, y - 0.95)

    # Central carbon
    draw_atom(ax, "C", x, y - 1.2)
    draw_bond(ax, x, y - 0.5, x, y - 1.2)

    # Double-bonded oxygen
    draw_atom(ax, "O", x + 0.7, y - 1.7)
    draw_bond(ax, x, y - 1.2, x + 0.7, y - 1.7)
    draw_bond(ax, x + 0.15, y - 1.2, x + 0.85, y - 1.7)  # parallel bond

    # Nitrogen
    draw_atom(ax, "N", x - 0.5, y - 1.6)
    draw_bond(ax, x, y - 1.2, x - 0.5, y - 1.6)

    # Two hydrogens on nitrogen
    draw_atom(ax, "H", x - 0.7, y - 2.0)
    draw_bond(ax, x - 0.5, y - 1.6, x - 0.7, y - 2.0)

    draw_atom(ax, "H", x - 0.3, y - 2.0)
    draw_bond(ax, x - 0.5, y - 1.6, x - 0.3, y - 2.0)

def draw_gln(ax, x, y):
    draw_atom(ax, "C", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)

    draw_atom(ax, "H", x - 0.3, y - 0.9)
    draw_bond(ax, x, y - 0.5, x - 0.3, y - 0.9)

    draw_atom(ax, "H", x + 0.3, y - 0.9)
    draw_bond(ax, x, y - 0.5, x + 0.3, y - 0.9)

    draw_atom(ax, "C", x, y - 1.2)
    draw_bond(ax, x, y - 0.5, x, y - 1.2)

    draw_atom(ax, "H", x - 0.3, y - 1.5)
    draw_bond(ax, x, y - 1.2, x - 0.3, y - 1.5)

    draw_atom(ax, "H", x + 0.3, y - 1.5)
    draw_bond(ax, x, y - 1.2, x + 0.3, y - 1.5)

    draw_atom(ax, "C", x, y - 1.9)
    draw_bond(ax, x, y - 1.2, x, y - 1.9)

    draw_atom(ax, "O", x + 0.5, y - 2.2)
    draw_bond(ax, x, y - 1.9, x + 0.5, y - 2.2)               # First bond line
    draw_bond(ax, x + 0.12, y - 1.9, x + 0.62, y - 2.2)       # Second (wider spaced) line

    draw_atom(ax, "N", x - 0.5, y - 2.2)
    draw_bond(ax, x, y - 1.9, x - 0.5, y - 2.2)

    draw_atom(ax, "H", x - 0.7, y - 2.6)
    draw_bond(ax, x - 0.5, y - 2.2, x - 0.7, y - 2.6)

    draw_atom(ax, "H", x - 0.3, y - 2.6)
    draw_bond(ax, x - 0.5, y - 2.2, x - 0.3, y - 2.6)

def draw_pro(ax, x, y):
    # Define the three ring carbon positions
    left = (x - 0.5, y - 0.7)
    bottom = (x, y - 1.3)
    right = (x + 0.5, y - 0.7)

    # Carbon atoms
    draw_atom(ax, "C", *left)
    draw_atom(ax, "C", *bottom)
    draw_atom(ax, "C", *right)

    # Ring bonds (forked connection at PRO base)
    draw_bond(ax, x - 0.2, y, *left)     # shifted left bond
    draw_bond(ax, *left, *bottom)
    draw_bond(ax, *bottom, *right)
    draw_bond(ax, *right, x + 0.2, y)    # shifted right bond

    # Hydrogens on each ring carbon (angled outward like your sketch)
    draw_atom(ax, "H", left[0] - 0.2, left[1])
    draw_bond(ax, *left, left[0] - 0.2, left[1])

    draw_atom(ax, "H", bottom[0], bottom[1] - 0.3)
    draw_bond(ax, *bottom, bottom[0], bottom[1] - 0.3)

    draw_atom(ax, "H", right[0] + 0.2, right[1])
    draw_bond(ax, *right, right[0] + 0.2, right[1])

def draw_asp(ax, x, y):
    # CH group
    draw_atom(ax, "C", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)

    draw_atom(ax, "H", x - 0.3, y - 0.9)
    draw_bond(ax, x, y - 0.5, x - 0.3, y - 0.9)

    draw_atom(ax, "H", x + 0.3, y - 0.9)
    draw_bond(ax, x, y - 0.5, x + 0.3, y - 0.9)

    # Central C
    draw_atom(ax, "C", x, y - 1.2)
    draw_bond(ax, x, y - 0.5, x, y - 1.2)

    # Double bonded O
    draw_atom(ax, "O", x - 0.5, y - 1.7)
    draw_bond(ax, x, y - 1.2, x - 0.5, y - 1.7)
    draw_bond(ax, x + 0.1, y - 1.2, x - 0.4, y - 1.7)  # wider parallel bond

    # Single bonded O⁻
    draw_atom(ax, "O⁻", x + 0.5, y - 1.7)
    draw_bond(ax, x, y - 1.2, x + 0.5, y - 1.7)

def draw_glu(ax, x, y):
    # CH group
    draw_atom(ax, "C", x, y - 0.5)
    draw_bond(ax, x, y, x, y - 0.5)

    draw_atom(ax, "H", x - 0.3, y - 0.9)
    draw_bond(ax, x, y - 0.5, x - 0.3, y - 0.9)

    draw_atom(ax, "H", x + 0.3, y - 0.9)
    draw_bond(ax, x, y - 0.5, x + 0.3, y - 0.9)

    # Middle CH2
    draw_atom(ax, "C", x, y - 1.2)
    draw_bond(ax, x, y - 0.5, x, y - 1.2)

    draw_atom(ax, "H", x - 0.3, y - 1.5)
    draw_bond(ax, x, y - 1.2, x - 0.3, y - 1.5)

    draw_atom(ax, "H", x + 0.3, y - 1.5)
    draw_bond(ax, x, y - 1.2, x + 0.3, y - 1.5)

    # Terminal C
    draw_atom(ax, "C", x, y - 1.9)
    draw_bond(ax, x, y - 1.2, x, y - 1.9)

    # Double bonded O
    draw_atom(ax, "O", x - 0.5, y - 2.4)
    draw_bond(ax, x, y - 1.9, x - 0.5, y - 2.4)
    draw_bond(ax, x + 0.1, y - 1.9, x - 0.4, y - 2.4)

    # Single bonded O⁻
    draw_atom(ax, "O⁻", x + 0.5, y - 2.4)
    draw_bond(ax, x, y - 1.9, x + 0.5, y - 2.4)


def draw_his(ax, x, y):
    r = 0.6          # Radius of imidazole ring
    dy = 1.3        # Vertical distance from CH2 to ring center
    offset = 0.08    # Double bond visual offset

    # CH2 below backbone node
    ch2_x, ch2_y = x, y - 0.5
    draw_atom(ax, "C", ch2_x, ch2_y)
    draw_bond(ax, x, y, ch2_x, ch2_y)

    # CH2 Hydrogens
    draw_atom(ax, "H", ch2_x - 0.3, ch2_y - 0.3)
    draw_bond(ax, ch2_x, ch2_y, ch2_x - 0.3, ch2_y - 0.3)
    draw_atom(ax, "H", ch2_x + 0.3, ch2_y - 0.3)
    draw_bond(ax, ch2_x, ch2_y, ch2_x + 0.3, ch2_y - 0.3)

    # Ring center (regular pentagon)
    ring_cx, ring_cy = ch2_x, ch2_y - dy

    # Generate regular upright pentagon starting from top (C1)
    ring_atoms = []
    for i in range(5):
        angle = math.pi / 2 + i * 2 * math.pi / 5  # 0° tilt
        px = ring_cx + r * math.cos(angle)
        py = ring_cy + r * math.sin(angle)
        ring_atoms.append((px, py))

    # Atom labels: C1, C2, N, C3, N⁺
    labels = ["C", "C", "N", "C", "N⁺"]
    ring = [(label, x, y) for label, (x, y) in zip(labels, ring_atoms)]

    # Draw ring atoms
    for label, px, py in ring:
        draw_atom(ax, label, px, py)

    # Draw ring bonds
    for i in range(5):
        x1, y1 = ring[i][1], ring[i][2]
        x2, y2 = ring[(i + 1) % 5][1], ring[(i + 1) % 5][2]
        draw_bond(ax, x1, y1, x2, y2)

    # Connect CH2 to C1 (index 0)
    draw_bond(ax, ch2_x, ch2_y, ring[0][1], ring[0][2])

    # === DOUBLE BONDS ===
    # C1 == C2 (index 0 & 1)
    x1, y1 = ring[0][1], ring[0][2]
    x2, y2 = ring[1][1], ring[1][2]
    dx, dy_ = x2 - x1, y2 - y1
    norm = math.hypot(dx, dy_)
    ox, oy = -dy_ / norm * offset, dx / norm * offset
    draw_bond(ax, x1 + ox, y1 + oy, x2 + ox, y2 + oy)

    # C3 == N⁺ (index 3 & 4)
    x3, y3 = ring[3][1], ring[3][2]
    x4, y4 = ring[4][1], ring[4][2]
    dx2, dy2 = x4 - x3, y4 - y3
    norm2 = math.hypot(dx2, dy2)
    ox2, oy2 = -dy2 / norm2 * offset, dx2 / norm2 * offset
    draw_bond(ax, x3 + ox2, y3 + oy2, x4 + ox2, y4 + oy2)

    # === Hydrogen on C2 (index 1) — left + down
    hx_c2 = ring[1][1] - 0.15
    hy_c2 = ring[1][2] - 0.3
    draw_atom(ax, "H", hx_c2, hy_c2)
    draw_bond(ax, ring[1][1], ring[1][2], hx_c2, hy_c2)

    # === Hydrogen on N (index 2) — left + down
    hx_n = ring[2][1] - 0.15
    hy_n = ring[2][2] - 0.3
    draw_atom(ax, "H", hx_n, hy_n)
    draw_bond(ax, ring[2][1], ring[2][2], hx_n, hy_n)

    # === Hydrogen on N⁺ (index 4) — right + down
    hx_nplus = ring[4][1] + 0.15
    hy_nplus = ring[4][2] - 0.3
    draw_atom(ax, "H", hx_nplus, hy_nplus)
    draw_bond(ax, ring[4][1], ring[4][2], hx_nplus, hy_nplus)


def draw_lys(ax, x, y):
    v = 0.7  # uniform vertical spacing

    # --- CH2 (from AA node) ---
    c1 = (x, y - 0.5)
    draw_atom(ax, "C", *c1)
    draw_bond(ax, x, y, *c1)

    draw_atom(ax, "H", x - 0.25, c1[1] - 0.3)
    draw_bond(ax, *c1, x - 0.25, c1[1] - 0.3)
    draw_atom(ax, "H", x + 0.25, c1[1] - 0.3)
    draw_bond(ax, *c1, x + 0.25, c1[1] - 0.3)

    # --- C2 ---
    c2 = (x, c1[1] - v)
    draw_atom(ax, "C", *c2)
    draw_bond(ax, *c1, *c2)
    draw_atom(ax, "H", x - 0.25, c2[1] - 0.3)
    draw_bond(ax, *c2, x - 0.25, c2[1] - 0.3)
    draw_atom(ax, "H", x + 0.25, c2[1] - 0.3)
    draw_bond(ax, *c2, x + 0.25, c2[1] - 0.3)

    # --- C3 ---
    c3 = (x, c2[1] - v)
    draw_atom(ax, "C", *c3)
    draw_bond(ax, *c2, *c3)
    draw_atom(ax, "H", x - 0.25, c3[1] - 0.3)
    draw_bond(ax, *c3, x - 0.25, c3[1] - 0.3)
    draw_atom(ax, "H", x + 0.25, c3[1] - 0.3)
    draw_bond(ax, *c3, x + 0.25, c3[1] - 0.3)

    # --- C4 ---
    c4 = (x, c3[1] - v)
    draw_atom(ax, "C", *c4)
    draw_bond(ax, *c3, *c4)
    draw_atom(ax, "H", x - 0.25, c4[1] - 0.3)
    draw_bond(ax, *c4, x - 0.25, c4[1] - 0.3)
    draw_atom(ax, "H", x + 0.25, c4[1] - 0.3)
    draw_bond(ax, *c4, x + 0.25, c4[1] - 0.3)

    # --- N⁺ ---
    n = (x, c4[1] - v)
    draw_atom(ax, "N⁺", *n)
    draw_bond(ax, *c4, *n)

    # NH₃⁺ Hydrogens (fanned out ~120°)
    angles = [math.radians(140), math.radians(80), math.radians(40)]
    for angle in angles:
        hx = n[0] + 0.5 * math.cos(angle)
        hy = n[1] - 0.5 * math.sin(angle)
        draw_atom(ax, "H", hx, hy)
        draw_bond(ax, *n, hx, hy)

import math

def draw_arg(ax, x, y):
    v = 0.7  # vertical step

    # CH2-1 (from node)
    c1 = (x, y - 0.5)
    draw_atom(ax, "C", *c1)
    draw_bond(ax, x, y, *c1)

    draw_atom(ax, "H", x - 0.25, c1[1] - 0.3)
    draw_bond(ax, *c1, x - 0.25, c1[1] - 0.3)
    draw_atom(ax, "H", x + 0.25, c1[1] - 0.3)
    draw_bond(ax, *c1, x + 0.25, c1[1] - 0.3)

    # CH2-2
    c2 = (x, c1[1] - v)
    draw_atom(ax, "C", *c2)
    draw_bond(ax, *c1, *c2)
    draw_atom(ax, "H", x - 0.25, c2[1] - 0.3)
    draw_bond(ax, *c2, x - 0.25, c2[1] - 0.3)
    draw_atom(ax, "H", x + 0.25, c2[1] - 0.3)
    draw_bond(ax, *c2, x + 0.25, c2[1] - 0.3)

    # CH2-3
    c3 = (x, c2[1] - v)
    draw_atom(ax, "C", *c3)
    draw_bond(ax, *c2, *c3)
    draw_atom(ax, "H", x + 0.25, c3[1] - 0.3)
    draw_bond(ax, *c3, x + 0.25, c3[1] - 0.3)

   # NH (middle N) – shift slightly right
    n_left = (x - 0.35, c3[1] - 0.35)
    draw_atom(ax, "N", *n_left)
    draw_bond(ax, *c3, *n_left)

    draw_atom(ax, "H", n_left[0] - 0.3, n_left[1])
    draw_bond(ax, *n_left, n_left[0] - 0.3, n_left[1])

    # Central guanidinium carbon – same vertical drop
    c_g = (x, n_left[1] - 0.4)
    draw_atom(ax, "C", *c_g)
    draw_bond(ax, *n_left, *c_g)

    # NH2 (left) – pushed down
    n2 = (x - 0.45, c_g[1] - 0.5)
    draw_atom(ax, "N", *n2)
    draw_bond(ax, *c_g, *n2)

    draw_atom(ax, "H", n2[0] - 0.2, n2[1] - 0.2)
    draw_bond(ax, *n2, n2[0] - 0.2, n2[1] - 0.2)

    draw_atom(ax, "H", n2[0] + 0.2, n2[1] - 0.2)
    draw_bond(ax, *n2, n2[0] + 0.2, n2[1] - 0.2)

    # NH2⁺ (right, double bonded) – also pushed down
    n3 = (x + 0.45, c_g[1] - 0.5)
    draw_atom(ax, "N⁺", *n3)
    draw_bond(ax, *c_g, *n3)

    # Double bond offset
    dx, dy = n3[0] - c_g[0], n3[1] - c_g[1]
    norm = math.hypot(dx, dy)
    ox, oy = -dy / norm * 0.08, dx / norm * 0.08
    draw_bond(ax, c_g[0] + ox, c_g[1] + oy, n3[0] + ox, n3[1] + oy)

    draw_atom(ax, "H", n3[0] + 0.2, n3[1] - 0.2)
    draw_bond(ax, *n3, n3[0] + 0.2, n3[1] - 0.2)

    draw_atom(ax, "H", n3[0] - 0.2, n3[1] - 0.2)
    draw_bond(ax, *n3, n3[0] - 0.2, n3[1] - 0.2)

def draw_trp(ax, x, y):
    # CH2 group below TRP node
    ch2 = (x, y - 0.5)
    draw_atom(ax, "C", *ch2)
    draw_bond(ax, x, y, *ch2)
    draw_atom(ax, "H", x - 0.25, ch2[1] - 0.3)
    draw_bond(ax, *ch2, x - 0.25, ch2[1] - 0.3)
    draw_atom(ax, "H", x + 0.25, ch2[1] - 0.3)
    draw_bond(ax, *ch2, x + 0.25, ch2[1] - 0.3)

    # Coordinates for indole ring (manual layout for visual accuracy)
    C1 = (x, y - 1.1)
    C2 = (x - 0.6, y - 1.3)
    C3 = (x - 0.7, y - 1.9)
    C4 = (x - 1.3, y - 2.1)
    C5 = (x - 1.4, y - 2.7)
    C6 = (x - 1.0, y - 3.2)
    C7 = (x - 0.4, y - 3.2)
    C8 = (x, y - 2.7)
    N  = (x + 0.5, y - 2.3)

    # Draw atoms
    for label, pos in [
        ("C", C1), ("C", C2), ("C", C3), ("C", C4),
        ("C", C5), ("C", C6), ("C", C7), ("C", C8), ("N", N)
    ]:
        draw_atom(ax, label, *pos)

    # Draw bonds (ring)
    draw_bond(ax, *ch2, *C1)
    draw_bond(ax, *C1, *C2)
    draw_bond(ax, *C2, *C3)
    draw_bond(ax, *C3, *C4)
    draw_bond(ax, *C4, *C5)
    draw_bond(ax, *C5, *C6)
    draw_bond(ax, *C6, *C7)
    draw_bond(ax, *C7, *C8)
    draw_bond(ax, *C8, *N)
    draw_bond(ax, *N, *C1)

    # Draw double bonds (aromatic rings)
    def dbl(a, b):
        dx, dy = b[0] - a[0], b[1] - a[1]
        norm = (dx**2 + dy**2)**0.5
        ox, oy = -dy / norm * 0.06, dx / norm * 0.06
        draw_bond(ax, a[0] + ox, a[1] + oy, b[0] + ox, b[1] + oy)

    dbl(C2, C3)
    dbl(C4, C5)
    dbl(C6, C7)
    dbl(C8, N)

    # Hydrogens on ring
    draw_atom(ax, "H", C4[0] - 0.2, C4[1] - 0.3)
    draw_bond(ax, *C4, C4[0] - 0.2, C4[1] - 0.3)

    draw_atom(ax, "H", C5[0] - 0.2, C5[1] - 0.3)
    draw_bond(ax, *C5, C5[0] - 0.2, C5[1] - 0.3)

    draw_atom(ax, "H", C6[0], C6[1] - 0.3)
    draw_bond(ax, *C6, C6[0], C6[1] - 0.3)

    draw_atom(ax, "H", N[0] + 0.25, N[1] - 0.2)
    draw_bond(ax, *N, N[0] + 0.25, N[1] - 0.2)



def draw_amino_acids(amino_acids):
    fig, ax = plt.subplots(figsize=(max(10, len(amino_acids) * 1.6), 6))
    ax.axis('off')
    ax.set_xlim(-1, len(amino_acids) * 3)
    ax.set_ylim(-3.5, 6)
    ax.set_aspect('equal')

    spacing = 3
    for i, aa in enumerate(amino_acids):
        x = i * spacing
        y_node = 4.5

        # Backbone connection (BEHIND everything else)
        if i > 0:
            ax.plot([x - spacing, x], [y_node, y_node], color='black', linewidth=1.5, zorder=0)

        # Node circle
        ax.add_patch(patches.Circle((x, y_node), 0.4,
                        facecolor=AA_COLORS.get(aa, 'lightgray'),
                        edgecolor='black', zorder=2))
        ax.text(x, y_node, aa if aa != "???" else "??", fontsize=10, ha='center', va='center', fontweight='bold', zorder=3)

        # R group bond (skip for PRO, it's drawn inside draw_pro)
        if aa != "PRO":
            ax.plot([x, x], [y_node - 0.4, y_node - 0.7], color='black', zorder=1)


        # Adjusted offset (R group higher)
        y_offset = y_node - 0.4

        if aa == "GLY":
            draw_gly(ax, x, y_offset)
        elif aa == "ALA":
            draw_ala(ax, x, y_offset)
        elif aa == "CYS":
            draw_cys(ax, x, y_offset)
        elif aa == "MET":
            draw_met(ax, x, y_offset)
        elif aa == "PHE":
            draw_phe(ax, x, y_offset)
        elif aa == "LEU":
            draw_leu(ax, x, y_offset)
        elif aa == "ILE":
            draw_ile(ax, x, y_offset)
        elif aa == "VAL":
            draw_val(ax, x, y_offset)
        elif aa == "SER":
            draw_ser(ax, x, y_offset)
        elif aa == "THR":
            draw_thr(ax, x, y_offset)
        elif aa == "ASN":
            draw_asn(ax, x, y_offset)
        elif aa == "GLN":
            draw_gln(ax, x, y_offset)
        elif aa == "TYR":
            draw_tyr(ax, x, y_offset)
        elif aa == "TRP":
            draw_trp(ax, x, y_offset)
        elif aa == "ASP":
            draw_asp(ax, x, y_offset)
        elif aa == "GLU":
            draw_glu(ax, x, y_offset)
        elif aa == "HIS":
            draw_his(ax, x, y_offset)
        elif aa == "LYS":
            draw_lys(ax, x, y_offset)
        elif aa == "ARG":
            draw_arg(ax, x, y_offset)
        elif aa == "PRO":
            draw_pro(ax, x, y_offset)
        elif aa == "STOP":
            pass  # Or draw_stop(ax, x, y_offset) if you want a symbol


    plt.title("Protein Structure Visualizer (Detailed Skeletal View)")
    plt.tight_layout()
    plt.show()

def validate_sequence(seq):
    seq = seq.upper().replace(" ", "")
    if not all(c in "ATCG" for c in seq):
        return False, "Invalid characters in sequence."
    if "ATG" not in seq:
        return False, "Missing START codon (ATG)."
    if not any(stop in seq for stop in ["TAA", "TAG", "TGA"]):
        return False, "Missing STOP codon."
    return True, seq

def split_into_codons(seq):
    return [seq[i:i+3] for i in range(0, len(seq), 3)]

def translate_codons(codons):
    result = []
    start = False
    for codon in codons:
        aa = CODON_TABLE.get(codon, "???")
        if aa == "MET":
            start = True
        if start:
            if aa == "STOP":
                break
            result.append(aa)
    return result

if __name__ == "__main__":
    dna_input = pyautogui.prompt("Enter DNA Sequence (e.g. ATGGCTTGTGGGTTTAA):")
    if dna_input:
        valid, cleaned = validate_sequence(dna_input)
        if not valid:
            pyautogui.alert("Error: " + cleaned)
        else:
            codons = split_into_codons(cleaned)
            aa_chain = translate_codons(codons)
            if not aa_chain:
                pyautogui.alert("No valid amino acids found.")
            else:
                draw_amino_acids(aa_chain)

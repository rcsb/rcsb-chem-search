# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause

"""
Periodic table of elements, without dependencies.
"""
from dataclasses import dataclass, field
from enum import Enum, auto, StrEnum, Flag
from functools import total_ordering
from typing import Self


class Category(StrEnum):
    ALKALI_METAL = auto()
    ALKALINE_EARTH_METAL = auto()
    LANTHANIDE = auto()
    ACTINIDE = auto()
    TRANSITION_METAL = auto()
    OTHER_METAL = auto()
    METALLOID = auto()
    OTHER_NONMETAL = auto()
    HALOGEN = auto()
    NOBLE_GAS = auto()


class Block(StrEnum):
    S = auto()
    P = auto()
    D = auto()
    F = auto()


@dataclass(frozen=True, slots=True)
class ElectronConfig:
    s: int
    p: int
    d: int
    f: int


class BondingType(Enum):
    metallic = auto()
    network_covalent = auto()
    molecular_covalent = auto()
    nobel = auto()
    unknown = auto()


@total_ordering
@dataclass(frozen=True, slots=True)
class Element:
    name: str
    symbol: str
    number: int
    type: Category
    period: int
    group: int
    block: Block
    valence: int
    oxidation_states: list[int]
    electron_config: ElectronConfig
    electronegativity: float
    bonding_type: BondingType
    weight_dalton: float
    atomic_radius_pm: float
    covalent_radius_pm: float
    van_der_waals_radius_pm: float

    def __lt__(self: Self, other: Self) -> bool:
        return self.atomic_number < other.atomic_number


class Elements(Enum):
    HYDROGEN = Element(
        name="Hydrogen",
        symbol="H",
        number=1,
        type=Category.OTHER_NONMETAL,
        period=1,
        group=1,
        block=Block.S,
        valence=1,
        oxidation_states=[-1, 1],
        electron_config=ElectronConfig(s=1, p=0, d=0, f=0),
        electronegativity=2.20,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=1.008,
        atomic_radius_pm=53,
        covalent_radius_pm=31,
        van_der_waals_radius_pm=120
    )

    HELIUM = Element(
        name="Helium",
        symbol="He",
        number=2,
        type=Category.NOBLE_GAS,
        period=1,
        group=18,
        block=Block.S,
        valence=0,
        oxidation_states=[],
        electron_config=ElectronConfig(s=2, p=0, d=0, f=0),
        electronegativity=0.0,
        bonding_type=BondingType.nobel,
        weight_dalton=4.002602,
        atomic_radius_pm=31,
        covalent_radius_pm=28,
        van_der_waals_radius_pm=140
    )

    LITHIUM = Element(
        name="Lithium",
        symbol="Li",
        number=3,
        type=Category.ALKALI_METAL,
        period=2,
        group=1,
        block=Block.S,
        valence=1,
        oxidation_states=[1],
        electron_config=ElectronConfig(s=2, p=1, d=0, f=0),
        electronegativity=0.98,
        bonding_type=BondingType.metallic,
        weight_dalton=6.94,
        atomic_radius_pm=167,
        covalent_radius_pm=128,
        van_der_waals_radius_pm=182
    )

    BERYLLIUM = Element(
        name="Beryllium",
        symbol="Be",
        number=4,
        type=Category.ALKALINE_EARTH_METAL,
        period=2,
        group=2,
        block=Block.S,
        valence=2,
        oxidation_states=[2],
        electron_config=ElectronConfig(s=2, p=2, d=0, f=0),
        electronegativity=1.57,
        bonding_type=BondingType.metallic,
        weight_dalton=9.0122,
        atomic_radius_pm=112,
        covalent_radius_pm=96,
        van_der_waals_radius_pm=153
    )

    BORON = Element(
        name="Boron",
        symbol="B",
        number=5,
        type=Category.METALLOID,
        period=2,
        group=13,
        block=Block.P,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=3, d=0, f=0),
        electronegativity=2.04,
        bonding_type=BondingType.network_covalent,
        weight_dalton=10.81,
        atomic_radius_pm=87,
        covalent_radius_pm=84,
        van_der_waals_radius_pm=192
    )

    CARBON = Element(
        name="Carbon",
        symbol="C",
        number=6,
        type=Category.OTHER_NONMETAL,
        period=2,
        group=14,
        block=Block.P,
        valence=4,
        oxidation_states=[-4, 2, 4],
        electron_config=ElectronConfig(s=2, p=4, d=0, f=0),
        electronegativity=2.55,
        bonding_type=BondingType.network_covalent,
        weight_dalton=12.011,
        atomic_radius_pm=67,
        covalent_radius_pm=77,
        van_der_waals_radius_pm=170
    )

    NITROGEN = Element(
        name="Nitrogen",
        symbol="N",
        number=7,
        type=Category.OTHER_NONMETAL,
        period=2,
        group=15,
        block=Block.P,
        valence=5,
        oxidation_states=[-3, 3, 5],
        electron_config=ElectronConfig(s=2, p=5, d=0, f=0),
        electronegativity=3.04,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=14.007,
        atomic_radius_pm=56,
        covalent_radius_pm=70,
        van_der_waals_radius_pm=155
    )

    OXYGEN = Element(
        name="Oxygen",
        symbol="O",
        number=8,
        type=Category.OTHER_NONMETAL,
        period=2,
        group=16,
        block=Block.P,
        valence=6,
        oxidation_states=[-2],
        electron_config=ElectronConfig(s=2, p=6, d=0, f=0),
        electronegativity=3.44,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=15.999,
        atomic_radius_pm=48,
        covalent_radius_pm=66,
        van_der_waals_radius_pm=152
    )

    FLUORINE = Element(
        name="Fluorine",
        symbol="F",
        number=9,
        type=Category.HALOGEN,
        period=2,
        group=17,
        block=Block.P,
        valence=7,
        oxidation_states=[-1],
        electron_config=ElectronConfig(s=2, p=7, d=0, f=0),
        electronegativity=3.98,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=18.998,
        atomic_radius_pm=42,
        covalent_radius_pm=64,
        van_der_waals_radius_pm=147
    )

    NEON = Element(
        name="Neon",
        symbol="Ne",
        number=10,
        type=Category.NOBLE_GAS,
        period=2,
        group=18,
        block=Block.P,
        valence=0,
        oxidation_states=[],
        electron_config=ElectronConfig(s=2, p=8, d=0, f=0),
        electronegativity=0.0,
        bonding_type=BondingType.nobel,
        weight_dalton=20.180,
        atomic_radius_pm=38,
        covalent_radius_pm=58,
        van_der_waals_radius_pm=154
    )

    SODIUM = Element(
        name="Sodium",
        symbol="Na",
        number=11,
        type=Category.ALKALI_METAL,
        period=3,
        group=1,
        block=Block.S,
        valence=1,
        oxidation_states=[1],
        electron_config=ElectronConfig(s=2, p=1, d=0, f=0),
        electronegativity=0.93,
        bonding_type=BondingType.metallic,
        weight_dalton=22.990,
        atomic_radius_pm=190,
        covalent_radius_pm=166,
        van_der_waals_radius_pm=227
    )

    MAGNESIUM = Element(
        name="Magnesium",
        symbol="Mg",
        number=12,
        type=Category.ALKALINE_EARTH_METAL,
        period=3,
        group=2,
        block=Block.S,
        valence=2,
        oxidation_states=[2],
        electron_config=ElectronConfig(s=2, p=2, d=0, f=0),
        electronegativity=1.31,
        bonding_type=BondingType.metallic,
        weight_dalton=24.305,
        atomic_radius_pm=145,
        covalent_radius_pm=141,
        van_der_waals_radius_pm=173
    )

    ALUMINIUM = Element(
        name="Aluminium",
        symbol="Al",
        number=13,
        type=Category.OTHER_METAL,
        period=3,
        group=13,
        block=Block.P,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=0, f=0),
        electronegativity=1.61,
        bonding_type=BondingType.metallic,
        weight_dalton=26.982,
        atomic_radius_pm=118,
        covalent_radius_pm=121,
        van_der_waals_radius_pm=184
    )

    SILICON = Element(
        name="Silicon",
        symbol="Si",
        number=14,
        type=Category.METALLOID,
        period=3,
        group=14,
        block=Block.P,
        valence=4,
        oxidation_states=[-4, 4],
        electron_config=ElectronConfig(s=2, p=6, d=0, f=0),
        electronegativity=1.90,
        bonding_type=BondingType.network_covalent,
        weight_dalton=28.085,
        atomic_radius_pm=111,
        covalent_radius_pm=111,
        van_der_waals_radius_pm=210
    )

    PHOSPHORUS = Element(
        name="Phosphorus",
        symbol="P",
        number=15,
        type=Category.OTHER_NONMETAL,
        period=3,
        group=15,
        block=Block.P,
        valence=5,
        oxidation_states=[-3, 3, 5],
        electron_config=ElectronConfig(s=2, p=6, d=0, f=0),
        electronegativity=2.19,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=30.974,
        atomic_radius_pm=98,
        covalent_radius_pm=107,
        van_der_waals_radius_pm=180
    )

    SULFUR = Element(
        name="Sulfur",
        symbol="S",
        number=16,
        type=Category.OTHER_NONMETAL,
        period=3,
        group=16,
        block=Block.P,
        valence=6,
        oxidation_states=[-2, 2, 4, 6],
        electron_config=ElectronConfig(s=2, p=6, d=0, f=0),
        electronegativity=2.58,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=32.06,
        atomic_radius_pm=88,
        covalent_radius_pm=105,
        van_der_waals_radius_pm=180
    )

    CHLORINE = Element(
        name="Chlorine",
        symbol="Cl",
        number=17,
        type=Category.HALOGEN,
        period=3,
        group=17,
        block=Block.P,
        valence=7,
        oxidation_states=[-1, 1, 3, 5, 7],
        electron_config=ElectronConfig(s=2, p=7, d=0, f=0),
        electronegativity=3.16,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=35.45,
        atomic_radius_pm=79,
        covalent_radius_pm=102,
        van_der_waals_radius_pm=175
    )

    ARGON = Element(
        name="Argon",
        symbol="Ar",
        number=18,
        type=Category.NOBLE_GAS,
        period=3,
        group=18,
        block=Block.P,
        valence=0,
        oxidation_states=[],
        electron_config=ElectronConfig(s=2, p=8, d=0, f=0),
        electronegativity=0.0,
        bonding_type=BondingType.nobel,
        weight_dalton=39.948,
        atomic_radius_pm=71,
        covalent_radius_pm=106,
        van_der_waals_radius_pm=188
    )

    POTASSIUM = Element(
        name="Potassium",
        symbol="K",
        number=19,
        type=Category.ALKALI_METAL,
        period=4,
        group=1,
        block=Block.S,
        valence=1,
        oxidation_states=[1],
        electron_config=ElectronConfig(s=2, p=6, d=0, f=0),
        electronegativity=0.82,
        bonding_type=BondingType.metallic,
        weight_dalton=39.098,
        atomic_radius_pm=243,
        covalent_radius_pm=203,
        van_der_waals_radius_pm=275
    )

    CALCIUM = Element(
        name="Calcium",
        symbol="Ca",
        number=20,
        type=Category.ALKALINE_EARTH_METAL,
        period=4,
        group=2,
        block=Block.S,
        valence=2,
        oxidation_states=[2],
        electron_config=ElectronConfig(s=2, p=6, d=0, f=0),
        electronegativity=1.00,
        bonding_type=BondingType.metallic,
        weight_dalton=40.078,
        atomic_radius_pm=194,
        covalent_radius_pm=176,
        van_der_waals_radius_pm=231
    )

    SCANDIUM = Element(
        name="Scandium",
        symbol="Sc",
        number=21,
        type=Category.TRANSITION_METAL,
        period=4,
        group=3,
        block=Block.D,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=0),
        electronegativity=1.36,
        bonding_type=BondingType.metallic,
        weight_dalton=44.956,
        atomic_radius_pm=184,
        covalent_radius_pm=170,
        van_der_waals_radius_pm=211
    )

    TITANIUM = Element(
        name="Titanium",
        symbol="Ti",
        number=22,
        type=Category.TRANSITION_METAL,
        period=4,
        group=4,
        block=Block.D,
        valence=4,
        oxidation_states=[2, 3, 4],
        electron_config=ElectronConfig(s=2, p=6, d=2, f=0),
        electronegativity=1.54,
        bonding_type=BondingType.metallic,
        weight_dalton=47.867,
        atomic_radius_pm=176,
        covalent_radius_pm=160,
        van_der_waals_radius_pm=210
    )

    VANADIUM = Element(
        name="Vanadium",
        symbol="V",
        number=23,
        type=Category.TRANSITION_METAL,
        period=4,
        group=5,
        block=Block.D,
        valence=5,
        oxidation_states=[2, 3, 4, 5],
        electron_config=ElectronConfig(s=2, p=6, d=3, f=0),
        electronegativity=1.63,
        bonding_type=BondingType.metallic,
        weight_dalton=50.942,
        atomic_radius_pm=171,
        covalent_radius_pm=153,
        van_der_waals_radius_pm=207
    )

    CHROMIUM = Element(
        name="Chromium",
        symbol="Cr",
        number=24,
        type=Category.TRANSITION_METAL,
        period=4,
        group=6,
        block=Block.D,
        valence=6,
        oxidation_states=[2, 3, 6],
        electron_config=ElectronConfig(s=2, p=6, d=5, f=0),
        electronegativity=1.66,
        bonding_type=BondingType.metallic,
        weight_dalton=51.996,
        atomic_radius_pm=166,
        covalent_radius_pm=139,
        van_der_waals_radius_pm=206
    )

    MANGANESE = Element(
        name="Manganese",
        symbol="Mn",
        number=25,
        type=Category.TRANSITION_METAL,
        period=4,
        group=7,
        block=Block.D,
        valence=7,
        oxidation_states=[2, 3, 4, 6, 7],
        electron_config=ElectronConfig(s=2, p=6, d=5, f=0),
        electronegativity=1.55,
        bonding_type=BondingType.metallic,
        weight_dalton=54.938,
        atomic_radius_pm=161,
        covalent_radius_pm=139,
        van_der_waals_radius_pm=205
    )

    IRON = Element(
        name="Iron",
        symbol="Fe",
        number=26,
        type=Category.TRANSITION_METAL,
        period=4,
        group=8,
        block=Block.D,
        valence=8,
        oxidation_states=[2, 3],
        electron_config=ElectronConfig(s=2, p=6, d=6, f=0),
        electronegativity=1.83,
        bonding_type=BondingType.metallic,
        weight_dalton=55.845,
        atomic_radius_pm=156,
        covalent_radius_pm=132,
        van_der_waals_radius_pm=204
    )

    COBALT = Element(
        name="Cobalt",
        symbol="Co",
        number=27,
        type=Category.TRANSITION_METAL,
        period=4,
        group=9,
        block=Block.D,
        valence=9,
        oxidation_states=[2, 3],
        electron_config=ElectronConfig(s=2, p=6, d=7, f=0),
        electronegativity=1.88,
        bonding_type=BondingType.metallic,
        weight_dalton=58.933,
        atomic_radius_pm=152,
        covalent_radius_pm=126,
        van_der_waals_radius_pm=204
    )

    NICKEL = Element(
        name="Nickel",
        symbol="Ni",
        number=28,
        type=Category.TRANSITION_METAL,
        period=4,
        group=10,
        block=Block.D,
        valence=10,
        oxidation_states=[2, 3],
        electron_config=ElectronConfig(s=2, p=6, d=8, f=0),
        electronegativity=1.91,
        bonding_type=BondingType.metallic,
        weight_dalton=58.693,
        atomic_radius_pm=149,
        covalent_radius_pm=124,
        van_der_waals_radius_pm=204
    )

    COPPER = Element(
        name="Copper",
        symbol="Cu",
        number=29,
        type=Category.TRANSITION_METAL,
        period=4,
        group=11,
        block=Block.D,
        valence=11,
        oxidation_states=[1, 2],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=0),
        electronegativity=1.90,
        bonding_type=BondingType.metallic,
        weight_dalton=63.546,
        atomic_radius_pm=145,
        covalent_radius_pm=132,
        van_der_waals_radius_pm=196
    )

    ZINC = Element(
        name="Zinc",
        symbol="Zn",
        number=30,
        type=Category.TRANSITION_METAL,
        period=4,
        group=12,
        block=Block.D,
        valence=12,
        oxidation_states=[2],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=0),
        electronegativity=1.65,
        bonding_type=BondingType.metallic,
        weight_dalton=65.38,
        atomic_radius_pm=142,
        covalent_radius_pm=122,
        van_der_waals_radius_pm=210
    )

    GALLIUM = Element(
        name="Gallium",
        symbol="Ga",
        number=31,
        type=Category.OTHER_METAL,
        period=4,
        group=13,
        block=Block.P,
        valence=3,
        oxidation_states=[1, 3],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=0),
        electronegativity=1.81,
        bonding_type=BondingType.metallic,
        weight_dalton=69.723,
        atomic_radius_pm=136,
        covalent_radius_pm=122,
        van_der_waals_radius_pm=187
    )

    GERMANIUM = Element(
        name="Germanium",
        symbol="Ge",
        number=32,
        type=Category.METALLOID,
        period=4,
        group=14,
        block=Block.P,
        valence=4,
        oxidation_states=[2, 4],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=0),
        electronegativity=2.01,
        bonding_type=BondingType.metallic,
        weight_dalton=72.630,
        atomic_radius_pm=125,
        covalent_radius_pm=122,
        van_der_waals_radius_pm=211
    )

    ARSENIC = Element(
        name="Arsenic",
        symbol="As",
        number=33,
        type=Category.METALLOID,
        period=4,
        group=15,
        block=Block.P,
        valence=5,
        oxidation_states=[-3, 3, 5],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=0),
        electronegativity=2.18,
        bonding_type=BondingType.network_covalent,
        weight_dalton=74.922,
        atomic_radius_pm=114,
        covalent_radius_pm=119,
        van_der_waals_radius_pm=185
    )

    SELENIUM = Element(
        name="Selenium",
        symbol="Se",
        number=34,
        type=Category.OTHER_NONMETAL,
        period=4,
        group=16,
        block=Block.P,
        valence=6,
        oxidation_states=[-2, 4, 6],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=0),
        electronegativity=2.55,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=78.971,
        atomic_radius_pm=103,
        covalent_radius_pm=120,
        van_der_waals_radius_pm=190
    )

    BROMINE = Element(
        name="Bromine",
        symbol="Br",
        number=35,
        type=Category.HALOGEN,
        period=4,
        group=17,
        block=Block.P,
        valence=7,
        oxidation_states=[-1, 1, 3, 5, 7],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=0),
        electronegativity=2.96,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=79.904,
        atomic_radius_pm=94,
        covalent_radius_pm=120,
        van_der_waals_radius_pm=185
    )

    KRYPTON = Element(
        name="Krypton",
        symbol="Kr",
        number=36,
        type=Category.NOBLE_GAS,
        period=4,
        group=18,
        block=Block.P,
        valence=0,
        oxidation_states=[0, 2],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=0),
        electronegativity=3.00,
        bonding_type=BondingType.nobel,
        weight_dalton=83.798,
        atomic_radius_pm=88,
        covalent_radius_pm=116,
        van_der_waals_radius_pm=202
    )

    RUBIDIUM = Element(
        name="Rubidium",
        symbol="Rb",
        number=37,
        type=Category.ALKALI_METAL,
        period=5,
        group=1,
        block=Block.S,
        valence=1,
        oxidation_states=[1],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=0),
        electronegativity=0.82,
        bonding_type=BondingType.metallic,
        weight_dalton=85.468,
        atomic_radius_pm=265,
        covalent_radius_pm=235,
        van_der_waals_radius_pm=303
    )

    STRONTIUM = Element(
        name="Strontium",
        symbol="Sr",
        number=38,
        type=Category.ALKALINE_EARTH_METAL,
        period=5,
        group=2,
        block=Block.S,
        valence=2,
        oxidation_states=[2],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=0),
        electronegativity=0.95,
        bonding_type=BondingType.metallic,
        weight_dalton=87.62,
        atomic_radius_pm=219,
        covalent_radius_pm=200,
        van_der_waals_radius_pm=249
    )

    YTTRIUM = Element(
        name="Yttrium",
        symbol="Y",
        number=39,
        type=Category.TRANSITION_METAL,
        period=5,
        group=3,
        block=Block.D,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=1),
        electronegativity=1.22,
        bonding_type=BondingType.metallic,
        weight_dalton=88.906,
        atomic_radius_pm=212,
        covalent_radius_pm=190,
        van_der_waals_radius_pm=212
    )

    ZIRCONIUM = Element(
        name="Zirconium",
        symbol="Zr",
        number=40,
        type=Category.TRANSITION_METAL,
        period=5,
        group=4,
        block=Block.D,
        valence=4,
        oxidation_states=[2, 3, 4],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=2),
        electronegativity=1.33,
        bonding_type=BondingType.metallic,
        weight_dalton=91.224,
        atomic_radius_pm=206,
        covalent_radius_pm=175,
        van_der_waals_radius_pm=206
    )

    NIOBIUM = Element(
        name="Niobium",
        symbol="Nb",
        number=41,
        type=Category.TRANSITION_METAL,
        period=5,
        group=5,
        block=Block.D,
        valence=5,
        oxidation_states=[2, 3, 5],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=3),
        electronegativity=1.60,
        bonding_type=BondingType.metallic,
        weight_dalton=92.906,
        atomic_radius_pm=198,
        covalent_radius_pm=164,
        van_der_waals_radius_pm=207
    )

    MOLYBDENUM = Element(
        name="Molybdenum",
        symbol="Mo",
        number=42,
        type=Category.TRANSITION_METAL,
        period=5,
        group=6,
        block=Block.D,
        valence=6,
        oxidation_states=[2, 3, 4, 5, 6],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=4),
        electronegativity=2.16,
        bonding_type=BondingType.metallic,
        weight_dalton=95.95,
        atomic_radius_pm=190,
        covalent_radius_pm=154,
        van_der_waals_radius_pm=209
    )

    TECHNETIUM = Element(
        name="Technetium",
        symbol="Tc",
        number=43,
        type=Category.TRANSITION_METAL,
        period=5,
        group=7,
        block=Block.D,
        valence=7,
        oxidation_states=[4, 7],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=5),
        electronegativity=1.90,
        bonding_type=BondingType.metallic,
        weight_dalton=98.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=183,
        covalent_radius_pm=147,
        van_der_waals_radius_pm=209
    )

    RUTHENIUM = Element(
        name="Ruthenium",
        symbol="Ru",
        number=44,
        type=Category.TRANSITION_METAL,
        period=5,
        group=8,
        block=Block.D,
        valence=8,
        oxidation_states=[3, 4, 6, 8],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=6),
        electronegativity=2.20,
        bonding_type=BondingType.metallic,
        weight_dalton=101.07,
        atomic_radius_pm=178,
        covalent_radius_pm=146,
        van_der_waals_radius_pm=207
    )

    RHODIUM = Element(
        name="Rhodium",
        symbol="Rh",
        number=45,
        type=Category.TRANSITION_METAL,
        period=5,
        group=9,
        block=Block.D,
        valence=9,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=7),
        electronegativity=2.28,
        bonding_type=BondingType.metallic,
        weight_dalton=102.91,
        atomic_radius_pm=173,
        covalent_radius_pm=142,
        van_der_waals_radius_pm=202
    )

    PALLADIUM = Element(
        name="Palladium",
        symbol="Pd",
        number=46,
        type=Category.TRANSITION_METAL,
        period=5,
        group=10,
        block=Block.D,
        valence=10,
        oxidation_states=[2, 4],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=8),
        electronegativity=2.20,
        bonding_type=BondingType.metallic,
        weight_dalton=106.42,
        atomic_radius_pm=169,
        covalent_radius_pm=139,
        van_der_waals_radius_pm=202
    )

    SILVER = Element(
        name="Silver",
        symbol="Ag",
        number=47,
        type=Category.TRANSITION_METAL,
        period=5,
        group=11,
        block=Block.D,
        valence=11,
        oxidation_states=[1],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=9),
        electronegativity=1.93,
        bonding_type=BondingType.metallic,
        weight_dalton=107.87,
        atomic_radius_pm=165,
        covalent_radius_pm=145,
        van_der_waals_radius_pm=172
    )

    CADMIUM = Element(
        name="Cadmium",
        symbol="Cd",
        number=48,
        type=Category.TRANSITION_METAL,
        period=5,
        group=12,
        block=Block.D,
        valence=12,
        oxidation_states=[2],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=10),
        electronegativity=1.69,
        bonding_type=BondingType.metallic,
        weight_dalton=112.41,
        atomic_radius_pm=161,
        covalent_radius_pm=144,
        van_der_waals_radius_pm=158
    )

    INDIUM = Element(
        name="Indium",
        symbol="In",
        number=49,
        type=Category.OTHER_METAL,
        period=5,
        group=13,
        block=Block.P,
        valence=3,
        oxidation_states=[1, 3],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=10),
        electronegativity=1.78,
        bonding_type=BondingType.metallic,
        weight_dalton=114.82,
        atomic_radius_pm=156,
        covalent_radius_pm=142,
        van_der_waals_radius_pm=193
    )

    TIN = Element(
        name="Tin",
        symbol="Sn",
        number=50,
        type=Category.OTHER_METAL,
        period=5,
        group=14,
        block=Block.P,
        valence=4,
        oxidation_states=[2, 4],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=10),
        electronegativity=1.96,
        bonding_type=BondingType.metallic,
        weight_dalton=118.71,
        atomic_radius_pm=145,
        covalent_radius_pm=139,
        van_der_waals_radius_pm=217
    )

    ANTIMONY = Element(
        name="Antimony",
        symbol="Sb",
        number=51,
        type=Category.METALLOID,
        period=5,
        group=15,
        block=Block.P,
        valence=5,
        oxidation_states=[-3, 3, 5],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=10),
        electronegativity=2.05,
        bonding_type=BondingType.network_covalent,
        weight_dalton=121.76,
        atomic_radius_pm=133,
        covalent_radius_pm=139,
        van_der_waals_radius_pm=206
    )

    TELLURIUM = Element(
        name="Tellurium",
        symbol="Te",
        number=52,
        type=Category.METALLOID,
        period=5,
        group=16,
        block=Block.P,
        valence=6,
        oxidation_states=[-2, 4, 6],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=10),
        electronegativity=2.10,
        bonding_type=BondingType.network_covalent,
        weight_dalton=127.60,
        atomic_radius_pm=123,
        covalent_radius_pm=138,
        van_der_waals_radius_pm=206
    )

    IODINE = Element(
        name="Iodine",
        symbol="I",
        number=53,
        type=Category.HALOGEN,
        period=5,
        group=17,
        block=Block.P,
        valence=7,
        oxidation_states=[-1, 1, 3, 5, 7],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=10),
        electronegativity=2.66,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=126.90,
        atomic_radius_pm=115,
        covalent_radius_pm=139,
        van_der_waals_radius_pm=198
    )

    XENON = Element(
        name="Xenon",
        symbol="Xe",
        number=54,
        type=Category.NOBLE_GAS,
        period=5,
        group=18,
        block=Block.P,
        valence=0,
        oxidation_states=[0, 2, 4, 6, 8],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=10),
        electronegativity=2.60,
        bonding_type=BondingType.nobel,
        weight_dalton=131.29,
        atomic_radius_pm=108,
        covalent_radius_pm=140,
        van_der_waals_radius_pm=216
    )

    CESIUM = Element(
        name="Cesium",
        symbol="Cs",
        number=55,
        type=Category.ALKALI_METAL,
        period=6,
        group=1,
        block=Block.S,
        valence=1,
        oxidation_states=[1],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=10),
        electronegativity=0.79,
        bonding_type=BondingType.metallic,
        weight_dalton=132.91,
        atomic_radius_pm=298,
        covalent_radius_pm=244,
        van_der_waals_radius_pm=343
    )

    BARIUM = Element(
        name="Barium",
        symbol="Ba",
        number=56,
        type=Category.ALKALINE_EARTH_METAL,
        period=6,
        group=2,
        block=Block.S,
        valence=2,
        oxidation_states=[2],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=10),
        electronegativity=0.89,
        bonding_type=BondingType.metallic,
        weight_dalton=137.33,
        atomic_radius_pm=253,
        covalent_radius_pm=215,
        van_der_waals_radius_pm=268
    )

    LANTHANUM = Element(
        name="Lanthanum",
        symbol="La",
        number=57,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=0),
        electronegativity=1.10,
        bonding_type=BondingType.metallic,
        weight_dalton=138.91,
        atomic_radius_pm=195,
        covalent_radius_pm=207,
        van_der_waals_radius_pm=240
    )

    CERIUM = Element(
        name="Cerium",
        symbol="Ce",
        number=58,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=4,
        oxidation_states=[3, 4],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=1),
        electronegativity=1.12,
        bonding_type=BondingType.metallic,
        weight_dalton=140.12,
        atomic_radius_pm=185,
        covalent_radius_pm=204,
        van_der_waals_radius_pm=235
    )

    PRASEODYMIUM = Element(
        name="Praseodymium",
        symbol="Pr",
        number=59,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=2),
        electronegativity=1.13,
        bonding_type=BondingType.metallic,
        weight_dalton=140.91,
        atomic_radius_pm=247,
        covalent_radius_pm=203,
        van_der_waals_radius_pm=239
    )

    NEODYMIUM = Element(
        name="Neodymium",
        symbol="Nd",
        number=60,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=3),
        electronegativity=1.14,
        bonding_type=BondingType.metallic,
        weight_dalton=144.24,
        atomic_radius_pm=206,
        covalent_radius_pm=201,
        van_der_waals_radius_pm=229
    )

    PROMETHIUM = Element(
        name="Promethium",
        symbol="Pm",
        number=61,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=4),
        electronegativity=1.13,
        bonding_type=BondingType.metallic,
        weight_dalton=145.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=205,
        covalent_radius_pm=199,
        van_der_waals_radius_pm=231
    )

    SAMARIUM = Element(
        name="Samarium",
        symbol="Sm",
        number=62,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[2, 3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=5),
        electronegativity=1.17,
        bonding_type=BondingType.metallic,
        weight_dalton=150.36,
        atomic_radius_pm=238,
        covalent_radius_pm=198,
        van_der_waals_radius_pm=229
    )

    EUROPIUM = Element(
        name="Europium",
        symbol="Eu",
        number=63,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[2, 3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=6),
        electronegativity=1.20,
        bonding_type=BondingType.metallic,
        weight_dalton=151.96,
        atomic_radius_pm=231,
        covalent_radius_pm=198,
        van_der_waals_radius_pm=233
    )

    GADOLINIUM = Element(
        name="Gadolinium",
        symbol="Gd",
        number=64,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=7),
        electronegativity=1.20,
        bonding_type=BondingType.metallic,
        weight_dalton=157.25,
        atomic_radius_pm=233,
        covalent_radius_pm=196,
        van_der_waals_radius_pm=237
    )

    TERBIUM = Element(
        name="Terbium",
        symbol="Tb",
        number=65,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=8),
        electronegativity=1.22,
        bonding_type=BondingType.metallic,
        weight_dalton=158.93,
        atomic_radius_pm=225,
        covalent_radius_pm=194,
        van_der_waals_radius_pm=221
    )

    DYSPROSIUM = Element(
        name="Dysprosium",
        symbol="Dy",
        number=66,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=9),
        electronegativity=1.23,
        bonding_type=BondingType.metallic,
        weight_dalton=162.50,
        atomic_radius_pm=228,
        covalent_radius_pm=192,
        van_der_waals_radius_pm=229
    )

    HOLMIUM = Element(
        name="Holmium",
        symbol="Ho",
        number=67,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=10),
        electronegativity=1.24,
        bonding_type=BondingType.metallic,
        weight_dalton=164.93,
        atomic_radius_pm=226,
        covalent_radius_pm=192,
        van_der_waals_radius_pm=218
    )

    ERBIUM = Element(
        name="Erbium",
        symbol="Er",
        number=68,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=11),
        electronegativity=1.24,
        bonding_type=BondingType.metallic,
        weight_dalton=167.26,
        atomic_radius_pm=226,
        covalent_radius_pm=189,
        van_der_waals_radius_pm=214
    )

    THULIUM = Element(
        name="Thulium",
        symbol="Tm",
        number=69,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[2, 3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=12),
        electronegativity=1.25,
        bonding_type=BondingType.metallic,
        weight_dalton=168.93,
        atomic_radius_pm=222,
        covalent_radius_pm=190,
        van_der_waals_radius_pm=215
    )

    YTTERBIUM = Element(
        name="Ytterbium",
        symbol="Yb",
        number=70,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=2,
        oxidation_states=[2, 3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=13),
        electronegativity=1.10,
        bonding_type=BondingType.metallic,
        weight_dalton=173.04,
        atomic_radius_pm=222,
        covalent_radius_pm=187,
        van_der_waals_radius_pm=219
    )

    LUTETIUM = Element(
        name="Lutetium",
        symbol="Lu",
        number=71,
        type=Category.LANTHANIDE,
        period=6,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=14),
        electronegativity=1.27,
        bonding_type=BondingType.metallic,
        weight_dalton=174.97,
        atomic_radius_pm=217,
        covalent_radius_pm=187,
        van_der_waals_radius_pm=217
    )

    HAFNIUM = Element(
        name="Hafnium",
        symbol="Hf",
        number=72,
        type=Category.TRANSITION_METAL,
        period=6,
        group=4,
        block=Block.D,
        valence=4,
        oxidation_states=[4],
        electron_config=ElectronConfig(s=2, p=6, d=2, f=14),
        electronegativity=1.30,
        bonding_type=BondingType.metallic,
        weight_dalton=178.49,
        atomic_radius_pm=208,
        covalent_radius_pm=175,
        van_der_waals_radius_pm=215
    )

    TANTALUM = Element(
        name="Tantalum",
        symbol="Ta",
        number=73,
        type=Category.TRANSITION_METAL,
        period=6,
        group=5,
        block=Block.D,
        valence=5,
        oxidation_states=[5],
        electron_config=ElectronConfig(s=2, p=6, d=3, f=14),
        electronegativity=1.50,
        bonding_type=BondingType.metallic,
        weight_dalton=180.95,
        atomic_radius_pm=200,
        covalent_radius_pm=170,
        van_der_waals_radius_pm=217
    )

    TUNGSTEN = Element(
        name="Tungsten",
        symbol="W",
        number=74,
        type=Category.TRANSITION_METAL,
        period=6,
        group=6,
        block=Block.D,
        valence=6,
        oxidation_states=[6],
        electron_config=ElectronConfig(s=2, p=6, d=4, f=14),
        electronegativity=2.36,
        bonding_type=BondingType.metallic,
        weight_dalton=183.84,
        atomic_radius_pm=193,
        covalent_radius_pm=162,
        van_der_waals_radius_pm=210
    )

    RHENIUM = Element(
        name="Rhenium",
        symbol="Re",
        number=75,
        type=Category.TRANSITION_METAL,
        period=6,
        group=7,
        block=Block.D,
        valence=7,
        oxidation_states=[2, 4, 6, 7],
        electron_config=ElectronConfig(s=2, p=6, d=5, f=14),
        electronegativity=1.90,
        bonding_type=BondingType.metallic,
        weight_dalton=186.21,
        atomic_radius_pm=188,
        covalent_radius_pm=151,
        van_der_waals_radius_pm=207
    )

    OSMIUM = Element(
        name="Osmium",
        symbol="Os",
        number=76,
        type=Category.TRANSITION_METAL,
        period=6,
        group=8,
        block=Block.D,
        valence=8,
        oxidation_states=[2, 3, 4, 6, 8],
        electron_config=ElectronConfig(s=2, p=6, d=6, f=14),
        electronegativity=2.20,
        bonding_type=BondingType.metallic,
        weight_dalton=190.23,
        atomic_radius_pm=185,
        covalent_radius_pm=144,
        van_der_waals_radius_pm=207
    )

    IRIDIUM = Element(
        name="Iridium",
        symbol="Ir",
        number=77,
        type=Category.TRANSITION_METAL,
        period=6,
        group=9,
        block=Block.D,
        valence=9,
        oxidation_states=[3, 4, 6],
        electron_config=ElectronConfig(s=2, p=6, d=7, f=14),
        electronegativity=2.20,
        bonding_type=BondingType.metallic,
        weight_dalton=192.22,
        atomic_radius_pm=180,
        covalent_radius_pm=141,
        van_der_waals_radius_pm=202
    )

    PLATINUM = Element(
        name="Platinum",
        symbol="Pt",
        number=78,
        type=Category.TRANSITION_METAL,
        period=6,
        group=10,
        block=Block.D,
        valence=10,
        oxidation_states=[2, 4],
        electron_config=ElectronConfig(s=2, p=6, d=9, f=14),
        electronegativity=2.28,
        bonding_type=BondingType.metallic,
        weight_dalton=195.08,
        atomic_radius_pm=177,
        covalent_radius_pm=136,
        van_der_waals_radius_pm=175
    )

    GOLD = Element(
        name="Gold",
        symbol="Au",
        number=79,
        type=Category.TRANSITION_METAL,
        period=6,
        group=11,
        block=Block.D,
        valence=11,
        oxidation_states=[1, 3],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=2.54,
        bonding_type=BondingType.metallic,
        weight_dalton=196.97,
        atomic_radius_pm=174,
        covalent_radius_pm=136,
        van_der_waals_radius_pm=166
    )

    MERCURY = Element(
        name="Mercury",
        symbol="Hg",
        number=80,
        type=Category.TRANSITION_METAL,
        period=6,
        group=12,
        block=Block.D,
        valence=12,
        oxidation_states=[1, 2],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=2.00,
        bonding_type=BondingType.metallic,
        weight_dalton=200.59,
        atomic_radius_pm=171,
        covalent_radius_pm=132,
        van_der_waals_radius_pm=155
    )

    THALLIUM = Element(
        name="Thallium",
        symbol="Tl",
        number=81,
        type=Category.OTHER_METAL,
        period=6,
        group=13,
        block=Block.P,
        valence=3,
        oxidation_states=[1, 3],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=1.62,
        bonding_type=BondingType.metallic,
        weight_dalton=204.38,
        atomic_radius_pm=156,
        covalent_radius_pm=145,
        van_der_waals_radius_pm=196
    )

    LEAD = Element(
        name="Lead",
        symbol="Pb",
        number=82,
        type=Category.OTHER_METAL,
        period=6,
        group=14,
        block=Block.P,
        valence=4,
        oxidation_states=[2, 4],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=2.33,
        bonding_type=BondingType.metallic,
        weight_dalton=207.2,
        atomic_radius_pm=154,
        covalent_radius_pm=146,
        van_der_waals_radius_pm=202
    )

    BISMUTH = Element(
        name="Bismuth",
        symbol="Bi",
        number=83,
        type=Category.OTHER_METAL,
        period=6,
        group=15,
        block=Block.P,
        valence=5,
        oxidation_states=[3, 5],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=2.02,
        bonding_type=BondingType.metallic,
        weight_dalton=208.98,
        atomic_radius_pm=143,
        covalent_radius_pm=148,
        van_der_waals_radius_pm=207
    )

    POLONIUM = Element(
        name="Polonium",
        symbol="Po",
        number=84,
        type=Category.METALLOID,
        period=6,
        group=16,
        block=Block.P,
        valence=6,
        oxidation_states=[2, 4, 6],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=2.00,
        bonding_type=BondingType.metallic,
        weight_dalton=209.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=135,
        covalent_radius_pm=140,
        van_der_waals_radius_pm=197
    )

    ASTATINE = Element(
        name="Astatine",
        symbol="At",
        number=85,
        type=Category.HALOGEN,
        period=6,
        group=17,
        block=Block.P,
        valence=7,
        oxidation_states=[-1, 1, 3, 5],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=2.20,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=210.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=127,
        covalent_radius_pm=150,
        van_der_waals_radius_pm=202
    )

    RADON = Element(
        name="Radon",
        symbol="Rn",
        number=86,
        type=Category.NOBLE_GAS,
        period=6,
        group=18,
        block=Block.P,
        valence=0,
        oxidation_states=[0],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=2.20,
        bonding_type=BondingType.nobel,
        weight_dalton=222.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=120,
        covalent_radius_pm=150,
        van_der_waals_radius_pm=220
    )

    FRANCIUM = Element(
        name="Francium",
        symbol="Fr",
        number=87,
        type=Category.ALKALI_METAL,
        period=7,
        group=1,
        block=Block.S,
        valence=1,
        oxidation_states=[1],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=0.7,
        bonding_type=BondingType.metallic,
        weight_dalton=223.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=348,
        covalent_radius_pm=260,
        van_der_waals_radius_pm=348
    )

    RADIUM = Element(
        name="Radium",
        symbol="Ra",
        number=88,
        type=Category.ALKALINE_EARTH_METAL,
        period=7,
        group=2,
        block=Block.S,
        valence=2,
        oxidation_states=[2],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=0.9,
        bonding_type=BondingType.metallic,
        weight_dalton=226.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=283,
        covalent_radius_pm=221,
        van_der_waals_radius_pm=283
    )

    ACTINIUM = Element(
        name="Actinium",
        symbol="Ac",
        number=89,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=3,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=0),
        electronegativity=1.1,
        bonding_type=BondingType.metallic,
        weight_dalton=227.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=295,
        covalent_radius_pm=215,
        van_der_waals_radius_pm=290
    )

    THORIUM = Element(
        name="Thorium",
        symbol="Th",
        number=90,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=4,
        oxidation_states=[4],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=1),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=232.04,
        atomic_radius_pm=285,
        covalent_radius_pm=206,
        van_der_waals_radius_pm=282
    )

    PROTACTINIUM = Element(
        name="Protactinium",
        symbol="Pa",
        number=91,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=5,
        oxidation_states=[4, 5],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=2),
        electronegativity=1.5,
        bonding_type=BondingType.metallic,
        weight_dalton=231.04,
        atomic_radius_pm=243,
        covalent_radius_pm=200,
        van_der_waals_radius_pm=243
    )

    URANIUM = Element(
        name="Uranium",
        symbol="U",
        number=92,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[3, 4, 5, 6],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=3),
        electronegativity=1.38,
        bonding_type=BondingType.metallic,
        weight_dalton=238.03,
        atomic_radius_pm=240,
        covalent_radius_pm=196,
        van_der_waals_radius_pm=240
    )

    NEPTUNIUM = Element(
        name="Neptunium",
        symbol="Np",
        number=93,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=7,
        oxidation_states=[3, 4, 5, 6, 7],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=4),
        electronegativity=1.36,
        bonding_type=BondingType.metallic,
        weight_dalton=237.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=221,
        covalent_radius_pm=196,
        van_der_waals_radius_pm=221
    )

    PLUTONIUM = Element(
        name="Plutonium",
        symbol="Pu",
        number=94,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[3, 4, 5, 6],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=6),
        electronegativity=1.28,
        bonding_type=BondingType.metallic,
        weight_dalton=244.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=243,
        covalent_radius_pm=187,
        van_der_waals_radius_pm=243
    )

    AMERICIUM = Element(
        name="Americium",
        symbol="Am",
        number=95,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[3, 4, 5, 6],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=7),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=243.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=244,
        covalent_radius_pm=180,
        van_der_waals_radius_pm=244
    )

    CURIUM = Element(
        name="Curium",
        symbol="Cm",
        number=96,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[3, 4],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=8),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=247.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=169,
        van_der_waals_radius_pm=245
    )

    BERKELIUM = Element(
        name="Berkelium",
        symbol="Bk",
        number=97,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[3, 4],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=9),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=247.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=244,
        covalent_radius_pm=166,
        van_der_waals_radius_pm=244
    )

    CALIFORNIUM = Element(
        name="Californium",
        symbol="Cf",
        number=98,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[3, 4],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=10),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=251.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=168,
        van_der_waals_radius_pm=245
    )

    EINSTEINIUM = Element(
        name="Einsteinium",
        symbol="Es",
        number=99,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=11),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=252.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=167,
        van_der_waals_radius_pm=245
    )

    FERMIUM = Element(
        name="Fermium",
        symbol="Fm",
        number=100,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=12),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=257.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=165,
        van_der_waals_radius_pm=245
    )

    MENDELEVIUM = Element(
        name="Mendelevium",
        symbol="Md",
        number=101,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[2, 3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=13),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=258.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=173,
        van_der_waals_radius_pm=245
    )

    NOBELIUM = Element(
        name="Nobelium",
        symbol="No",
        number=102,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[2, 3],
        electron_config=ElectronConfig(s=2, p=6, d=1, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=259.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=172,
        van_der_waals_radius_pm=245
    )

    LAWRENCIUM = Element(
        name="Lawrencium",
        symbol="Lr",
        number=103,
        type=Category.ACTINIDE,
        period=7,
        group=3,
        block=Block.F,
        valence=6,
        oxidation_states=[3],
        electron_config=ElectronConfig(s=2, p=6, d=2, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=262.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=172,
        van_der_waals_radius_pm=245
    )

    RUTHERFORDIUM = Element(
        name="Rutherfordium",
        symbol="Rf",
        number=104,
        type=Category.TRANSITION_METAL,
        period=7,
        group=4,
        block=Block.D,
        valence=4,
        oxidation_states=[4],
        electron_config=ElectronConfig(s=2, p=6, d=2, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=267.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=157,
        van_der_waals_radius_pm=245
    )

    DUBNIUM = Element(
        name="Dubnium",
        symbol="Db",
        number=105,
        type=Category.TRANSITION_METAL,
        period=7,
        group=5,
        block=Block.D,
        valence=5,
        oxidation_states=[5],
        electron_config=ElectronConfig(s=2, p=6, d=3, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=270.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=149,
        van_der_waals_radius_pm=245
    )

    SEABORGIUM = Element(
        name="Seaborgium",
        symbol="Sg",
        number=106,
        type=Category.TRANSITION_METAL,
        period=7,
        group=6,
        block=Block.D,
        valence=6,
        oxidation_states=[6],
        electron_config=ElectronConfig(s=2, p=6, d=4, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=271.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=143,
        van_der_waals_radius_pm=245
    )

    BOHRIUM = Element(
        name="Bohrium",
        symbol="Bh",
        number=107,
        type=Category.TRANSITION_METAL,
        period=7,
        group=7,
        block=Block.D,
        valence=7,
        oxidation_states=[7],
        electron_config=ElectronConfig(s=2, p=6, d=5, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=270.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=141,
        van_der_waals_radius_pm=245
    )

    HASSIUM = Element(
        name="Hassium",
        symbol="Hs",
        number=108,
        type=Category.TRANSITION_METAL,
        period=7,
        group=8,
        block=Block.D,
        valence=8,
        oxidation_states=[8],
        electron_config=ElectronConfig(s=2, p=6, d=6, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=277.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=134,
        van_der_waals_radius_pm=245
    )

    MEITNERIUM = Element(
        name="Meitnerium",
        symbol="Mt",
        number=109,
        type=Category.TRANSITION_METAL,
        period=7,
        group=9,
        block=Block.D,
        valence=9,
        oxidation_states=[9],
        electron_config=ElectronConfig(s=2, p=6, d=7, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=278.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=129,
        van_der_waals_radius_pm=245
    )

    DARMSTADTIUM = Element(
        name="Darmstadtium",
        symbol="Ds",
        number=110,
        type=Category.TRANSITION_METAL,
        period=7,
        group=10,
        block=Block.D,
        valence=10,
        oxidation_states=[10],
        electron_config=ElectronConfig(s=2, p=6, d=8, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=281.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=128,
        van_der_waals_radius_pm=245
    )

    ROENTGENIUM = Element(
        name="Roentgenium",
        symbol="Rg",
        number=111,
        type=Category.TRANSITION_METAL,
        period=7,
        group=11,
        block=Block.D,
        valence=11,
        oxidation_states=[11],
        electron_config=ElectronConfig(s=2, p=6, d=9, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=282.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=121,
        van_der_waals_radius_pm=245
    )

    COPERNICIUM = Element(
        name="Copernicium",
        symbol="Cn",
        number=112,
        type=Category.TRANSITION_METAL,
        period=7,
        group=12,
        block=Block.D,
        valence=12,
        oxidation_states=[12],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=285.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=122,
        van_der_waals_radius_pm=245
    )

    NIHONIUM = Element(
        name="Nihonium",
        symbol="Nh",
        number=113,
        type=Category.OTHER_METAL,
        period=7,
        group=13,
        block=Block.P,
        valence=3,
        oxidation_states=[1, 3],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=286.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=136,
        van_der_waals_radius_pm=245
    )

    FLEROVIUM = Element(
        name="Flerovium",
        symbol="Fl",
        number=114,
        type=Category.OTHER_METAL,
        period=7,
        group=14,
        block=Block.P,
        valence=4,
        oxidation_states=[2, 4],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=289.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=143,
        van_der_waals_radius_pm=245
    )

    MOSCOVIUM = Element(
        name="Moscovium",
        symbol="Mc",
        number=115,
        type=Category.OTHER_METAL,
        period=7,
        group=15,
        block=Block.P,
        valence=5,
        oxidation_states=[1, 3, 5],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=290.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=143,
        van_der_waals_radius_pm=245
    )

    LIVERMORIUM = Element(
        name="Livermorium",
        symbol="Lv",
        number=116,
        type=Category.OTHER_METAL,
        period=7,
        group=16,
        block=Block.P,
        valence=6,
        oxidation_states=[2, 4],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.metallic,
        weight_dalton=293.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=156,
        van_der_waals_radius_pm=245
    )

    TENNESSINE = Element(
        name="Tennessine",
        symbol="Ts",
        number=117,
        type=Category.HALOGEN,
        period=7,
        group=17,
        block=Block.P,
        valence=7,
        oxidation_states=[1, 3, 5, 7],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.molecular_covalent,
        weight_dalton=294.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=156,
        van_der_waals_radius_pm=245
    )

    OGANESSON = Element(
        name="Oganesson",
        symbol="Og",
        number=118,
        type=Category.NOBLE_GAS,
        period=7,
        group=18,
        block=Block.P,
        valence=0,
        oxidation_states=[0],
        electron_config=ElectronConfig(s=2, p=6, d=10, f=14),
        electronegativity=1.3,
        bonding_type=BondingType.nobel,
        weight_dalton=294.0,  # Radioactive, most stable isotope's mass
        atomic_radius_pm=245,
        covalent_radius_pm=156,
        van_der_waals_radius_pm=245
    )

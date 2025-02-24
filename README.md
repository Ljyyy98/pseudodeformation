# Pseudorepresentations and Deformation Rings

This repository was created by **Jinyue Luo** to provide computational tools supporting the paper  
["Pseudorepresentations not arising from genuine representations"](https://arxiv.org/abs/2310.16953).

## Overview

In this paper, we investigate the conditions under which a pseudorepresentation of a finite group \( G \) over a commutative ring \( A \) may fail to be realized by an actual group representation. This study is crucial for understanding the limitations of lifting pseudorepresentations to genuine representations, especially in the context of deformation theory and number theory.

## Contents

The repository includes the following scripts:

- **Magma Script (`pseudodeformation_ring.m`)**  
  Computes the pseudodeformation ring over \( \mathbb{F}_2 \) for a specified finite group using the Magma computational algebra system.

- **SageMath Script (`rcoeff_pseudodeformation.sage`)**  
  Evaluates specific elements within the framed deformation ring over the integer ring, leveraging SageMath for symbolic computation.

## Prerequisites

To utilize these scripts, ensure you have the following software installed:

- **[Magma](https://magma.maths.usyd.edu.au/)** – Required for running `pseudodeformation_ring.m`.  
- **[SageMath](https://www.sagemath.org/)** – Required for executing `rcoeff_pseudodeformation.sage`.


## Notes

- The **Magma script** computes the pseudodeformation ring over \( \mathbb{F}_2 \), using the group's GAP ID.
- The **SageMath script** evaluates specific elements within the framed deformation ring over the integers.

For more details, refer to the full paper:  
["Pseudorepresentations not arising from genuine representations"](https://arxiv.org/abs/2310.16953).


# BitBIRCH is an open-source clustering module based on iSIM
#
# Please, cite the BitBIRCH paper: https://www.biorxiv.org/content/10.1101/2024.08.10.607459v1
#
# BitBIRCH is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# BitBIRCH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# BitBIRCH authors (PYTHON): Ramon Alain Miranda Quintana <ramirandaq@gmail.com>, <quintana@chem.ufl.edu>
#                            Vicky (Vic) Jung <jungvicky@ufl.edu>
#                            Kenneth Lopez Perez <klopezperez@chem.ufl.edu>
#                            Kate Huddleston <kdavis2@chem.ufl.edu>
#
# BitBIRCH License: LGPL-3.0 https://www.gnu.org/licenses/lgpl-3.0.en.html#license-text

import bitbirch.bitbirch as bb


##Functions##
#add docstrings to functions
def radius(fps, threshold = 0.65, branching_factor=50):
    """Runnign radius birch; classic birch"""
    bb.set_merge('radius')
    brc = bb.BitBirch(branching_factor=branching_factor, threshold=threshold)
    brc.fit(fps)
    return brc

def diameter(fps, threshold = 0.65):
    bb.set_merge('diameter')
    brc = bb.BitBirch(branching_factor=50, threshold=threshold)
    brc.fit(fps)
    return brc

def diameter_prune(fps, threshold = 0.65):
    bb.set_merge('diameter')
    brc = bb.BitBirch(branching_factor=50, threshold=threshold)
    brc.fit(fps, singly=False)
    brc.prune(fps)
    return brc

def diameter_prune_tolerance(fps, tol=0.05, threshold = 0.65):
    bb.set_merge('diameter')
    brc = bb.BitBirch(branching_factor=50, threshold=threshold)
    brc.fit(fps, singly=False)
    bb.set_merge('tolerance', tolerance=tol)
    brc.prune(fps)
    return brc

def diameter_prune_reassign(fps, threshold = 0.65):
    bb.set_merge('diameter')
    brc = bb.BitBirch(branching_factor=50, threshold=threshold)  
    brc.fit(fps, singly=False)
    brc.prune(fps)
    brc.reassign(fps)
    return brc

def diameter_prune_tolerance_reassign(fps, tol=0.05, threshold = 0.65):
    bb.set_merge('diameter')
    brc = bb.BitBirch(branching_factor=50, threshold=threshold)
    brc.fit(fps, singly=False)
    bb.set_merge('tolerance', tolerance=tol)
    brc.prune(fps)
    brc.reassign(fps)
    return brc

def BFs_reclustering(fps, init_threshold = 0.65, second_threshold = 0.7, second_tolerance = 0.0):
    # Do the diameter clustering
    bb.set_merge('diameter')
    brc = bb.BitBirch(branching_factor=50, threshold=init_threshold)
    brc.fit(fps, singly=True)

    # Extract the BFs for the second clustering
    BFs = brc.prepare_BFs(fps)

    # Do the second clustering
    bb.set_merge('tolerance', tolerance=second_tolerance)
    brc = bb.BitBirch(branching_factor=50, threshold=second_threshold)
    brc.fit_BFs(BFs) # Note that we fit the BFs, not the fps

    return brc
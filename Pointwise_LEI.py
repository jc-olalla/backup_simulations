"""Pointwise glyph input file generator"""
import subprocess
import k_parametric_LEI
import numpy as np


def mesh_generation_pointiwse(glyph_script_path, mesh_file_path, Re, points):
    # Pointwise normal mesh extrusion
    initStep = k_parametric_LEI.wall_height(Re)   # Initial mesh layer height [m]
    print(initStep)
    growthFactor = 1.1                  # Growth factor
    dimTotal = 575                      # Total airfoil surface dimension
    n_TE = 20                           # TE dimension
    n_LE = int(dimTotal*0.3) #150   0.31  0.27perc                     # LE dimension
    n_upper = int((dimTotal-n_TE-n_LE)/2)   # upper surface dimension
    n_lower = n_upper                       # lower surface dimension
    numLayers = 201                    # Number of mesh layers (away from the wall)

    # Hyperbolic smoothing parameters
    explSmooth = 8.0                   #0.5 Explicit smoothing: Default is 0.5; between 0 and 10
    implSmooth = 2 * explSmooth        #1 Implicit smoothing: Default is 1.0; between 0 and infinity but double explicit smoothing
    kbSmooth = 5.0                     #0 Kinsey Barth: Default is 0.0; greater than 3 if mesh front includes severe concavities
    volSmooth = 0.5                    #0.5 Volume smoothing: Default is 0.5; between 0 and 1

    # Write the GLF file
    with open(glyph_script_path, "w") as f:
        f.write("package require PWI_Glyph 7.22.2\n")

        # Create each point with an identifier
        for i, point in enumerate(points):
            f.write(f"  set _DB({i+1}) [pw::Point create]\n")
            f.write(f"  $_DB({i+1}) setPoint {{{point[0]} {point[1]} {point[2]}}}\n")

        # Start creating the spline
        f.write("set _TMP(mode_1) [pw::Application begin Create]\n")
        f.write("  set _TMP(PW_1) [pw::SegmentSpline create]\n")

        # link identifiers to a spline
        for i, point in enumerate(points):
            f.write(f"  $_TMP(PW_1) addPoint [list 0 0 $_DB({i+1})]\n")
        f.write(f"\n")

        # Setup airfoil spline
        f.write("$_TMP(PW_1) setSlope Akima\n")
        f.write("set _CN(1) [pw::Connector create]\n")
        f.write("$_CN(1) addSegment $_TMP(PW_1)\n")
        f.write("$_CN(1) calculateDimension\n")
        f.write("unset _TMP(PW_1)\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n\n")

        # Point LE
        f.write("lappend _TMP(split_params) [$_CN(1) getParameter -closest [pw::Application getXYZ [$_CN(1) closestControlPoint [list 0 0 $_DB(20)]]]]\n")
        f.write("set _TMP(PW_1) [$_CN(1) split $_TMP(split_params)]\n")
        f.write("unset _TMP(PW_1)\n")
        f.write("unset _TMP(split_params)\n\n")

        # point fillet
        f.write("set _CN(2) [pw::GridEntity getByName con-1-split-2]\n")
        f.write("lappend _TMP(split_params) [$_CN(2) getParameter -closest [pw::Application getXYZ [$_CN(2) closestControlPoint [list 0 0 $_DB(370)]]]]\n")
        f.write("set _TMP(PW_1) [$_CN(2) split $_TMP(split_params)]\n")
        f.write("unset _TMP(PW_1)\n")
        f.write("unset _TMP(split_params)\n\n")

        # join
        f.write("set _CN(3) [pw::GridEntity getByName con-1-split-1]\n")
        f.write("set _CN(4) [pw::GridEntity getByName con-1-split-2-split-2]\n")
        f.write("set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(3) $_CN(4)]]\n")
        f.write("unset _TMP(ignored)\n")
        f.write("unset _TMP(PW_1)\n\n")

        # TE
        f.write("set _CN(5) [pw::GridEntity getByName con-1-split-2-split-1]\n")
        f.write("lappend _TMP(split_params) [$_CN(5) getParameter -closest [pw::Application getXYZ [$_CN(5) closestControlPoint [list 0 0 $_DB(179)]]]]\n")
        f.write("lappend _TMP(split_params) [$_CN(5) getParameter -closest [pw::Application getXYZ [$_CN(5) closestControlPoint [list 0 0 $_DB(208)]]]]\n")
        f.write("set _TMP(PW_1) [$_CN(5) split $_TMP(split_params)]\n")
        f.write("unset _TMP(PW_1)\n")
        f.write("unset _TMP(split_params)\n\n")


        # Set dimension of uniform distribution of TE
        f.write("set _CN(1) [pw::GridEntity getByName con-1-split-2-split-1-split-2]\n")
        f.write("set _TMP(mode_1) [pw::Application begin Dimension]\n")
        f.write("set _TMP(PW_1) [pw::Collection create]\n")
        f.write("$_TMP(PW_1) set [list $_CN(1)]\n")  # Add $_CN(3) to collection
        f.write(f"$_TMP(PW_1) do setDimension -resetDistribution {n_TE}\n")  # Set dimension distribution
        f.write("$_TMP(PW_1) delete\n")  # Delete collection after use
        f.write("unset _TMP(PW_1)\n")  # Clean up
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n\n")

        # Calculate spacing TE
        f.write(f"set uniform_spacing_TE [pwu::Vector3 length [pwu::Vector3 subtract [$_CN(1) getXYZ {2}] [$_CN(1) getXYZ {1}]]]\n")
        f.write("puts \"Spacing TE: $uniform_spacing_TE\"\n\n")


        # Set dimension of uniform distribution of LE
        f.write("set _CN(2) [pw::GridEntity getByName con-1-split-1]\n")
        f.write("set _TMP(mode_1) [pw::Application begin Dimension]\n")
        f.write("set _TMP(PW_1) [pw::Collection create]\n")
        f.write("$_TMP(PW_1) set [list $_CN(2)]\n")  # Add $_CN(3) to collection
        f.write(f"$_TMP(PW_1) do setDimension -resetDistribution {n_LE}\n")  # Set dimension distribution
        f.write("$_TMP(PW_1) delete\n")  # Delete collection after use
        f.write("unset _TMP(PW_1)\n")  # Clean up
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n\n")

        # Calculate spacing LE
        f.write(f"set uniform_spacing_LE [pwu::Vector3 length [pwu::Vector3 subtract [$_CN(2) getXYZ {2}] [$_CN(2) getXYZ {1}]]]\n")
        f.write("puts \"Spacing LE: $uniform_spacing_LE\"\n\n")

        # Distribution Top surface
        f.write("set _CN(3) [pw::GridEntity getByName con-1-split-2-split-1-split-1]\n")
        f.write("set _TMP(mode_1) [pw::Application begin Modify [list $_CN(3)]]\n")
        f.write("$_CN(3) replaceDistribution 1 [pw::DistributionTanh create]\n")
        f.write("[$_CN(3) getDistribution 1] setBeginSpacing $uniform_spacing_LE\n")
        f.write("[$_CN(3) getDistribution 1] setEndSpacing $uniform_spacing_TE\n")  # Set end spacing
        f.write("$_CN(3) setSubConnectorDimensionFromDistribution 1\n")  # Apply the distribution
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n\n")

        # Dimension Top surface
        f.write("set _TMP(mode_2) [pw::Application begin Dimension]\n")
        f.write("set _TMP(PW_1) [pw::Collection create]\n")
        f.write("$_TMP(PW_1) set [list $_CN(3)]\n")  # Add connector to collection
        f.write(f"$_TMP(PW_1) do setDimension -resetDistribution {n_upper}\n")  # Set dimension distribution
        f.write("$_TMP(PW_1) delete\n")  # Clean up collection
        f.write("unset _TMP(PW_1)\n")
        f.write("$_TMP(mode_2) end\n")
        f.write("unset _TMP(mode_2)\n\n")


        # Distribution lower surface
        f.write("set _CN(4) [pw::GridEntity getByName con-1-split-2-split-1-split-3]\n")
        f.write("set _TMP(mode_1) [pw::Application begin Modify [list $_CN(4)]]\n")
        f.write("$_CN(4) replaceDistribution 1 [pw::DistributionTanh create]\n")
        f.write("[$_CN(4) getDistribution 1] setBeginSpacing $uniform_spacing_TE\n")
        f.write("[$_CN(4) getDistribution 1] setEndSpacing $uniform_spacing_LE\n")  # Set end spacing
        f.write("$_CN(4) setSubConnectorDimensionFromDistribution 1\n")  # Apply the distribution
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n\n")

        # Dimension lower surface
        f.write("set _TMP(mode_2) [pw::Application begin Dimension]\n")
        f.write("set _TMP(PW_1) [pw::Collection create]\n")
        f.write("$_TMP(PW_1) set [list $_CN(4)]\n")  # Add connector to collection
        f.write(f"$_TMP(PW_1) do setDimension -resetDistribution {n_lower}\n")  # Set dimension distribution
        f.write("$_TMP(PW_1) delete\n")  # Clean up collection
        f.write("unset _TMP(PW_1)\n")
        f.write("$_TMP(mode_2) end\n")
        f.write("unset _TMP(mode_2)\n\n")


        # Normal Extrusion
        f.write("set _TMP(mode_1) [pw::Application begin Create]\n")
        f.write("set _TMP(PW_1) [pw::Edge createFromConnectors [list $_CN(1) $_CN(2) $_CN(3) $_CN(4)]]\n")
        f.write("set _TMP(edge_1) [lindex $_TMP(PW_1) 0]\n")
        f.write("unset _TMP(PW_1)\n")
        f.write("set _DM(1) [pw::DomainStructured create]\n")
        f.write("$_DM(1) addEdge $_TMP(edge_1)\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n")
        f.write("set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_DM(1)]]\n")
        f.write("$_TMP(mode_1) setKeepFailingStep false\n")
        f.write("$_DM(1) setExtrusionSolverAttribute SpacingMode Algebraic\n")
        f.write(f"$_DM(1) setExtrusionSolverAttribute NormalInitialStepSize {initStep}\n")
        f.write(f"$_DM(1) setExtrusionSolverAttribute SpacingGrowthFactor {growthFactor}\n")
        f.write("$_DM(1) setExtrusionSolverAttribute NormalMarchingVector {0 0 1}\n")
        f.write("$_DM(1) setExtrusionSolverAttribute StopAtPositiveSkewJacobian true\n")
        f.write("$_DM(1) setExtrusionSolverAttribute StopAtZeroJacobian true\n")
        f.write("$_DM(1) setExtrusionSolverAttribute StopAtNegativeJacobian true\n")
        f.write("$_DM(1) setExtrusionSolverAttribute StopAtNegativeSkewJacobian true\n")
        f.write("$_DM(1) setExtrusionSolverAttribute Mode NormalHyperbolic\n")
        f.write(f"$_DM(1) setExtrusionSolverAttribute NormalExplicitSmoothing {explSmooth}\n")
        f.write(f"$_DM(1) setExtrusionSolverAttribute NormalImplicitSmoothing {implSmooth}\n")
        f.write(f"$_DM(1) setExtrusionSolverAttribute NormalKinseyBarthSmoothing {kbSmooth}\n")
        f.write(f"$_DM(1) setExtrusionSolverAttribute NormalVolumeSmoothing {volSmooth}\n")
        f.write(f"$_TMP(mode_1) run {numLayers}\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n")
        f.write("unset _TMP(edge_1)\n\n")


        # Create 3D block by 1 translation step in z-direction
        f.write("set _TMP(mode_1) [pw::Application begin Create]\n")
        f.write("set _DM(1) [pw::GridEntity getByName dom-1]\n")
        f.write("set _TMP(PW_1) [pw::FaceStructured createFromDomains [list $_DM(1)]]\n")
        f.write("set _TMP(face_1) [lindex $_TMP(PW_1) 0]\n")
        f.write("unset _TMP(PW_1)\n")
        f.write("set _BL(1) [pw::BlockStructured create]\n")
        f.write("$_BL(1) addFace $_TMP(face_1)\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n\n")

        f.write("set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_BL(1)]]\n")
        f.write("$_TMP(mode_1) setKeepFailingStep true\n")
        f.write("$_BL(1) setExtrusionSolverAttribute Mode Translate\n")
        f.write("$_BL(1) setExtrusionSolverAttribute TranslateDirection {0 0 1}\n")
        f.write("$_BL(1) setExtrusionSolverAttribute TranslateDistance 1\n")
        f.write("$_TMP(mode_1) run 1\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n")
        f.write("unset _TMP(face_1)\n\n")


        # Boundary Conditions
        f.write("pw::Application setCAESolver OpenFOAM 3\n\n")


        # Farfield
        f.write("set _TMP(PW_1) [pw::BoundaryCondition create]\n")
        f.write("set _TMP(PW_1) [pw::BoundaryCondition getByName bc-2]\n")
        f.write("set _DM(2) [pw::GridEntity getByName dom-7]\n")
        f.write("$_TMP(PW_1) apply [list [list $_BL(1) $_DM(2)]]\n")
        f.write("$_TMP(PW_1) setName farfield\n")
        f.write("$_TMP(PW_1) setPhysicalType -usage CAE patch\n")
        f.write("unset _TMP(PW_1)\n\n")


        # Sides
        f.write("set _TMP(PW_2) [pw::BoundaryCondition create]\n")
        f.write("set _TMP(PW_2) [pw::BoundaryCondition getByName bc-3]\n")
        f.write("set _DM(3) [pw::GridEntity getByName dom-1]\n")
        f.write("set _DM(4) [pw::GridEntity getByName dom-9]\n")
        f.write("$_TMP(PW_2) apply [list [list $_BL(1) $_DM(3)] [list $_BL(1) $_DM(4)]]\n")
        f.write("$_TMP(PW_2) setName side\n")
        f.write("$_TMP(PW_2) setPhysicalType -usage CAE empty\n")
        f.write("unset _TMP(PW_2)\n\n")


        # Airfoil
        f.write("set _TMP(PW_3) [pw::BoundaryCondition create]\n")
        f.write("set _TMP(PW_3) [pw::BoundaryCondition getByName bc-4]\n")
        f.write("set _DM(5) [pw::GridEntity getByName dom-2]\n")
        f.write("set _DM(6) [pw::GridEntity getByName dom-3]\n")
        f.write("set _DM(7) [pw::GridEntity getByName dom-4]\n")
        f.write("set _DM(8) [pw::GridEntity getByName dom-5]\n")
        f.write("$_TMP(PW_3) apply [list [list $_BL(1) $_DM(5)] [list $_BL(1) $_DM(6)] [list $_BL(1) $_DM(7)] [list $_BL(1) $_DM(8)]]\n")
        f.write("$_TMP(PW_3) setName airfoil\n")
        f.write("$_TMP(PW_3) setPhysicalType -usage CAE wall\n")
        f.write("unset _TMP(PW_3)\n\n")


        # Export to CAE
        f.write("set _TMP(mode_1) [pw::Application begin CaeExport]\n")
        f.write("$_TMP(mode_1) addAllEntities\n")
        f.write(f"$_TMP(mode_1) initialize -strict -type CAE {mesh_file_path}\n")
        f.write("$_TMP(mode_1) verify\n")
        f.write("$_TMP(mode_1) write\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n\n")

        # Home screen
        f.write("pw::Display resetView -Z\n")

        # Exit application
        f.write("pw::Application exit\n")
    return ()
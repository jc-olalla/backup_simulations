"""Pointwise glyph input file generator"""
import subprocess
from toolbox import k_parametric_LEI
import numpy as np


class GlfScript():
    def __init__():
        pass

def initialize_glyph(glyph_script_path):
    # Write the GLF file
    with open(glyph_script_path, "w") as f:
        f.write("package require PWI_Glyph 8.24.1\n")

def create_database(glyph_script_path, df_points):
    df_points.sort_values(by='db', inplace=True)
    with open(glyph_script_path, "a") as f:
        for index, row in df_points.iterrows():
            # Create each point with an identifier
            #print(f"Creating _DB({row['db']})")
            f.write(f"set _DB({row['db']}) [pw::Point create]\n")
            f.write(f"  $_DB({row['db']}) setPoint {{{row['x']} {row['y']} 0.0}}\n")

def get_db_ids(df_points, *pieces):
    dbs = []
    for piece in pieces:
        df_piece_i = df_points.copy()[df_points['piece_id'] == piece]
        dbs.append(df_piece_i['db'].to_list())

    db_ids = sum(dbs, [])
    return db_ids


def create_connector_with_db_ids(glyph_script_path, db_ids, connector_id, slope_type='Akima'):
    """
    segment_type = 'SegmentLinear'
    creates connector
    connector_id (int)
    """
    with open(glyph_script_path, "a") as f:
        # Start creating the spline
        f.write("\n"*2)
        f.write("set _TMP(mode_1) [pw::Application begin Create]\n")
        f.write(f"set _TMP(PW_1) [pw::SegmentSpline create]\n")

        # link identifiers to a spline
        for point_id in db_ids:
            f.write(f"  $_TMP(PW_1) addPoint [list 0 0 $_DB({point_id})]\n")
        f.write(f"\n")

        # Setup airfoil spline
        f.write(f"$_TMP(PW_1) setSlope {slope_type}\n")
        f.write(f"set _CN({connector_id}) [pw::Connector create]\n")
        f.write(f"$_CN({connector_id}) addSegment $_TMP(PW_1)\n")
        f.write(f"$_CN({connector_id}) calculateDimension\n")
        f.write("unset _TMP(PW_1)\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n\n")

def add_dimension_to_connector(glyph_script_path, connector_id, ref_connector='con-1', dimension=50):
    # connector_id could be any number (I think)
    with open(glyph_script_path, "a") as f:
        f.write(f"set _CN({connector_id}) [pw::GridEntity getByName {ref_connector}]\n")
        f.write("set _TMP(mode_2) [pw::Application begin Dimension]\n")
        f.write("set _TMP(PW_1) [pw::Collection create]\n")
        f.write(f"$_TMP(PW_1) set [list $_CN({connector_id})]\n")
        f.write(f"$_TMP(PW_1) do setDimension -resetDistribution {dimension}\n")
        f.write("$_TMP(PW_1) delete\n")
        f.write("unset _TMP(PW_1)\n")
        f.write("$_TMP(mode_2) end\n")
        f.write("unset _TMP(mode_2)\n")

def join_connectors(glyph_script_path, *connector_ids):
    with open(glyph_script_path, "a") as f:
        f.write("\n"*2)
        con_str = '[list'
        for i, con_id in enumerate(connector_ids):
            con_str += f' $_CN({i+1})'
            f.write(f"set _CN({i+1}) [pw::GridEntity getByName {con_id}]\n")
        con_str += ']'

    with open(glyph_script_path, "a") as f:
        f.write(f"set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution {con_str}]\n")
        #f.write(f"set _CN({new_id}) [lindex $_TMP(PW_1) 0]\n")
        f.write("unset _TMP(ignored)\n")
        f.write("unset _TMP(PW_1)\n\n")

        # rename this joined connector set _CN(100)

#def create_domain_from_connectors(glyph_script_path, connector_id, ref_connector='con-4'):
#    with open(glyph_script_path, "a") as f:
#        f.write("\n"*2)
#        f.write("set _TMP(mode_1) [pw::Application begin Create]\n")
#        f.write(f"set _CN({connector_id}) [pw::GridEntity getByName {ref_connector}]\n")
#        f.write(f"set _TMP(PW_1) [pw::Edge createFromConnectors [list $_CN({connector_id})]]\n")
#        f.write("set _TMP(edge_1) [lindex $_TMP(PW_1) 0]\n")
#        f.write("unset _TMP(PW_1)\n")
#        f.write("set _DM(1) [pw::DomainStructured create]\n")
#        f.write("$_DM(1) addEdge $_TMP(edge_1)\n")
#        f.write("$_TMP(mode_1) end\n")
#        f.write("unset _TMP(mode_1)\n")

def create_domain_from_connectors(glyph_script_path, *connector_ids):
    with open(glyph_script_path, "a") as f:
        f.write("\n"*2)
        f.write("set _TMP(mode_1) [pw::Application begin Create]\n")
        con_str = ''
        for i, con_id in enumerate(connector_ids):
            con_str += f' $_CN({i+1})'
            f.write(f"set _CN({i+1}) [pw::GridEntity getByName {con_id}]\n")

        f.write(f"set _TMP(PW_1) [pw::Edge createFromConnectors [list {con_str}]]\n")
        f.write("set _TMP(edge_1) [lindex $_TMP(PW_1) 0]\n")
        f.write("unset _TMP(PW_1)\n")
        f.write("set _DM(1) [pw::DomainStructured create]\n")
        f.write("$_DM(1) addEdge $_TMP(edge_1)\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n")


def create_c_domain_from_connectors(glyph_script_path, connectors):
    with open(glyph_script_path, "a") as f:
        f.write("set _TMP(mode_1) [pw::Application begin Create]\n")
        unique_connectors = []
        mapper = {}
        for i, con_i in enumerate(connectors):
            if con_i not in unique_connectors:
                f.write(f"set _CN({i+1}) [pw::GridEntity getByName {con_i}]\n")
                unique_connectors.append(con_i)
                mapper[f'{con_i}'] = f'_CN({i+1})'

        f.write("set _TMP(edge_1) [pw::Edge create]\n")

        for i, con_i in enumerate(connectors):
            CN = mapper[f'{con_i}']
            f.write(f"$_TMP(edge_1) addConnector ${CN}\n")

        f.write("set _DM(1) [pw::DomainStructured create]\n")
        f.write("$_DM(1) addEdge $_TMP(edge_1)\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n")


def hyperbolic_extrusion(glyph_script_path, Re, growthFactor=1.1, explSmooth=0.5, implSmooth=1.0, kbSmooth=0.0, volSmooth=0.5, numLayers=201):
    initStep = k_parametric_LEI.wall_height(Re)   # Initial mesh layer height [m]
    print(f'Initial extrusion step: {initStep}')
    #initStep = 0.0021
    with open(glyph_script_path, "a") as f:
        # Normal Extrusion
        f.write("\n"*2)
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
        f.write("unset _TMP(edge_1)")

def hyperbolic_extrusion_case2(glyph_script_path, Re, growthFactor=1.1, explSmooth=0.5, implSmooth=1.0, kbSmooth=0.0, volSmooth=0.5, numLayers=201):
    #initStep = k_parametric_LEI.wall_height(Re)   # Initial mesh layer height [m]
    initStep = 0.001  # maybe 0.0001  or change smoothing
    print(f'Initial extrusion step: {initStep}')
    with open(glyph_script_path, "a") as f:
        # Normal Extrusion
        f.write("\n"*2)
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
        f.write("unset _TMP(edge_1)")

def double_hyperbolic_extrusion(glyph_script_path, Re, growthFactor=1.1, explSmooth=0.5, implSmooth=1.0, kbSmooth=0.0, volSmooth1=0.1, volSmooth2=0.5, numLayers1=10, numLayers2=100):
    initStep = k_parametric_LEI.wall_height(Re)   # Initial mesh layer height [m]
    with open(glyph_script_path, "a") as f:
        # Normal Extrusion
        f.write("\n"*2)
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
        f.write(f"$_DM(1) setExtrusionSolverAttribute NormalVolumeSmoothing {volSmooth1}\n")
        f.write(f"$_TMP(mode_1) run {numLayers1}\n")
        f.write(f"$_DM(1) setExtrusionSolverAttribute NormalVolumeSmoothing {volSmooth2}\n")
        f.write(f"$_TMP(mode_1) run {numLayers2}\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n")
        f.write("unset _TMP(edge_1)")



def extrusion_to_3d(glyph_script_path, *domain_ids):
    with open(glyph_script_path, "a") as f:
        # Create 3D block by 1 translation step in z-direction
        f.write("\n"*2)
        f.write("set _TMP(mode_1) [pw::Application begin Create]\n")
        dom_str = ''
        bl_str = ''
        for i, dom_i in enumerate(domain_ids):
            f.write(f"set _DM({i+1}) [pw::GridEntity getByName {dom_i}]\n")
            dom_str += f' $_DM({i+1})'
            bl_str += f' $_BL({i+1})'
        f.write(f"set _TMP(PW_1) [pw::FaceStructured createFromDomains [list{dom_str}]]\n")

        for i, dom_i in enumerate(domain_ids):
            f.write(f"set _TMP(face_{i+1}) [lindex $_TMP(PW_1) {i}]\n")
        f.write("unset _TMP(PW_1)\n")

        for i, dom_i in enumerate(domain_ids):
            f.write(f"set _BL({i+1}) [pw::BlockStructured create]\n")
            f.write(f"$_BL({i+1}) addFace $_TMP(face_{i+1})\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n\n")

        f.write(f"set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list{bl_str}]]\n")

        f.write("$_TMP(mode_1) setKeepFailingStep true\n")
        f.write("$_BL(1) setExtrusionSolverAttribute Mode Translate\n")
        f.write("$_BL(1) setExtrusionSolverAttribute TranslateDirection {0 0 1}\n")
        f.write("$_BL(1) setExtrusionSolverAttribute TranslateDistance 1\n")
        f.write("$_TMP(mode_1) run 1\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n")
        f.write("unset _TMP(face_1)\n\n")

def assign_boundary_conditions(glyph_script_path):
    # TODO: generalize for multiple airfoil components
    # Boundary Conditions
    with open(glyph_script_path, "a") as f:
        f.write("\n"*2)
        f.write("pw::Application setCAESolver OpenFOAM 3\n\n")

        # Farfield
        f.write("set _TMP(PW_1) [pw::BoundaryCondition create]\n")
        f.write("set _TMP(PW_1) [pw::BoundaryCondition getByName bc-2]\n")
        f.write("set _DM(2) [pw::GridEntity getByName dom-4]\n")
        f.write("$_TMP(PW_1) apply [list [list $_BL(1) $_DM(2)]]\n")
        f.write("$_TMP(PW_1) setName farfield\n")
        f.write("$_TMP(PW_1) setPhysicalType -usage CAE patch\n")
        f.write("unset _TMP(PW_1)\n\n")


        # Sides
        f.write("set _TMP(PW_2) [pw::BoundaryCondition create]\n")
        f.write("set _TMP(PW_2) [pw::BoundaryCondition getByName bc-3]\n")
        f.write("set _DM(3) [pw::GridEntity getByName dom-1]\n")
        f.write("set _DM(4) [pw::GridEntity getByName dom-6]\n")
        f.write("$_TMP(PW_2) apply [list [list $_BL(1) $_DM(3)] [list $_BL(1) $_DM(4)]]\n")
        f.write("$_TMP(PW_2) setName side\n")
        f.write("$_TMP(PW_2) setPhysicalType -usage CAE empty\n")
        f.write("unset _TMP(PW_2)\n\n")

        # Airfoil
        f.write("set _TMP(PW_3) [pw::BoundaryCondition create]\n")
        f.write("set _TMP(PW_3) [pw::BoundaryCondition getByName bc-4]\n")
        f.write("set _DM(5) [pw::GridEntity getByName dom-2]\n")
        f.write("$_TMP(PW_3) apply [list [list $_BL(1) $_DM(5)]]\n")
        f.write("$_TMP(PW_3) setName airfoil\n")
        f.write("$_TMP(PW_3) setPhysicalType -usage CAE wall\n")
        f.write("unset _TMP(PW_3)\n\n")

    #dom-1 mesh face 1
    #dom-2 airfoil surface
    #dom-3 2d domain used for creating the airfoil >> not part of the mesh
    #dom-4 airfoil-outer rings surface
    #dom-6 mesh face 2

def assign_boundary_conditions_case2(glyph_script_path):
    # TODO: generalize for multiple airfoil components
    # Boundary Conditions
    with open(glyph_script_path, "a") as f:
        f.write("\n"*2)
        f.write("pw::Application setCAESolver OpenFOAM 3\n\n")

        # Create farfield bc
        f.write("set _TMP(PW_1) [pw::BoundaryCondition create]\n")
        f.write("pw::Application markUndoLevel {Create BC}\n")
        f.write("unset _TMP(PW_1)\n")

        f.write("set _TMP(PW_1) [pw::BoundaryCondition getByName bc-2]\n")
        f.write("$_TMP(PW_1) setName farfield\n")
        f.write("pw::Application markUndoLevel {Name BC}\n")

        f.write("set _BL(1) [pw::GridEntity getByName blk-2]\n")  # outer mesh block
        f.write("set _DM(2) [pw::GridEntity getByName dom-17]\n")
        f.write("$_TMP(PW_1) apply [list [list $_BL(1) $_DM(2)]]\n")
        f.write("pw::Application markUndoLevel {Set BC}\n")

        f.write("$_TMP(PW_1) setPhysicalType -usage CAE patch\n")
        f.write("pw::Application markUndoLevel {Change BC type}\n")
        f.write("unset _TMP(PW_1)\n\n")

        # Create side bc
        f.write("set _TMP(PW_2) [pw::BoundaryCondition create]\n")
        f.write("pw::Application markUndoLevel {Create BC}\n")
        f.write("unset _TMP(PW_2)\n")

        f.write("set _TMP(PW_2) [pw::BoundaryCondition getByName bc-3]\n")
        f.write("$_TMP(PW_2) setName side\n")
        f.write("pw::Application markUndoLevel {Name BC}\n")

        f.write("set _BL(1) [pw::GridEntity getByName blk-2]\n")  # outer mesh block
        f.write("set _BL(2) [pw::GridEntity getByName blk-1]\n")  # disk mesh block
        f.write("set _DM(3) [pw::GridEntity getByName dom-1]\n")
        f.write("set _DM(4) [pw::GridEntity getByName dom-2]\n")
        f.write("set _DM(5) [pw::GridEntity getByName dom-9]\n")
        f.write("set _DM(6) [pw::GridEntity getByName dom-19]\n")
        f.write("$_TMP(PW_2) apply [list [list $_BL(2) $_DM(3)] [list $_BL(1) $_DM(4)] [list $_BL(2) $_DM(5)] [list $_BL(1) $_DM(6)]]\n")
        f.write("pw::Application markUndoLevel {Set BC}\n")

        f.write("$_TMP(PW_2) setPhysicalType -usage CAE empty\n")
        f.write("pw::Application markUndoLevel {Change BC type}\n")
        f.write("unset _TMP(PW_2)\n\n")

    #dom-1 air rotor mesh 1 
    #dom-2 o-grid mesh 1
    #dom-9 air rotor mesh 2
    #dom-17 farfield
    #dom-19 o-grid mesh 2

def export_mesh(glyph_script_path, mesh_file_path):
    with open(glyph_script_path, "a") as f:
        f.write("\n"*2)
        # Export to CAE
        f.write("set _BL(1) [pw::GridEntity getByName blk-1]")
        f.write("set _TMP(mode_1) [pw::Application begin CaeExport [pw::Entity sort [list $_BL(1)]]]\n")
        f.write(f"$_TMP(mode_1) initialize -strict -type CAE {{mesh_file_path}}\n")
        f.write("$_TMP(mode_1) verify\n")
        f.write("$_TMP(mode_1) write\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)")

def redistribute_connector_points(glyph_script_path, con_left='con-1', con_mid='con-2', con_right='con-3', matching_point_con_left='end', matching_point_con_right='beginning'):
    # TODO: check if renaming the connectors would affect the overal workflow of the script
    with open(glyph_script_path, "a") as f:
        # "Left" connector con_left
        f.write("\n")
        f.write(f"set _CN(1000) [pw::GridEntity getByName {con_left}]\n")
        if matching_point_con_left == 'end':
            f.write("set numPoints_left [$_CN(1000) getPointCount]\n")
            f.write("set pt_last_left     [$_CN(1000) getXYZ [expr {$numPoints_left - 1}]]\n")
            f.write("set pt_2nd_last_left [$_CN(1000) getXYZ [expr {$numPoints_left - 2}]]\n")
            f.write("set delta_ini [pwu::Vector3 length [pwu::Vector3 subtract $pt_last_left $pt_2nd_last_left]]\n")
        elif matching_point_con_left == 'beginning':
            f.write("set delta_ini [pwu::Vector3 length [pwu::Vector3 subtract [$_CN(1000) getXYZ 2] [$_CN(1000) getXYZ 1]]]\n")
        else:
            f.write(f"set delta_ini {matching_point_con_left}\n")  # hack to allow custom distance


        # "Right" connector con_right
        f.write("\n")
        f.write(f"set _CN(3000) [pw::GridEntity getByName {con_right}]\n")
        if matching_point_con_right == 'beginning':
            f.write("set delta_end [pwu::Vector3 length [pwu::Vector3 subtract [$_CN(3000) getXYZ 2] [$_CN(3000) getXYZ 1]]]\n")
        elif matching_point_con_right == 'end':
            f.write("numPoints_right [$_CN(3000) getPointCount]\n")
            f.write("set pt_last_right     [$_CN(3000) getXYZ [expr {$numPoints_right - 1}]]\n")
            f.write("set pt_2nd_last_right [$_CN(3000) getXYZ [expr {$numPoints_right - 2}]]\n")
            f.write("set delta_end [pwu::Vector3 length [pwu::Vector3 subtract $pt_last_right $pt_2nd_last_right]]\n")
            f.write(f"set delta_end {matching_point_con_left}\n")  # hack to allow custom distance


        # Redistribute connector points con_mid
        f.write("\n")
        f.write(f"set _CN(2000) [pw::GridEntity getByName {con_mid}]\n")
        f.write(f"set _TMP(mode_1) [pw::Application begin Modify [list $_CN(2000)]]\n")
        f.write(f"$_CN(2000) replaceDistribution 1 [pw::DistributionTanh create]\n")
        f.write(f"[$_CN(2000) getDistribution 1] setBeginSpacing $delta_ini\n")
        f.write(f"[$_CN(2000) getDistribution 1] setEndSpacing $delta_end\n")
        f.write(f"$_CN(2000) setSubConnectorDimensionFromDistribution 1\n")
        f.write(f"$_TMP(mode_1) end\n")
        f.write(f"unset _TMP(mode_1)\n")

def assemble_connectors(glyph_script_path, *connector_ids):
    with open(glyph_script_path, "a") as f:
        f.write("\n"*2)
        f.write(f"set _TMP(mode_1) [pw::Application begin Create]\n")
        con_str = ''
        for i, con_id in enumerate(connector_ids):
            con_str += f' $_CN({i+1})'
            f.write(f"set _CN({i+1}) [pw::GridEntity getByName {con_id}]\n")

        f.write(f"set _TMP(PW_1) [pw::DomainStructured createFromConnectors -reject _TMP(unusedCons) [list {con_str}]]\n")
        f.write(f"unset _TMP(unusedCons)\n")
        f.write("unset _TMP(PW_1)\n")
        f.write("$_TMP(mode_1) end\n")
        f.write("unset _TMP(mode_1)\n\n")

def assign_cell_zones(glyph_script_path, xmin, ymin, xmax, ymax, zone_name="actuatorDiskZone"):
    with open(glyph_script_path, "a") as f:
        f.write("\n"*2)
        f.write(f"set xmin {xmin}\n")
        f.write(f"set xmax {xmax}\n")
        f.write(f"set ymin {ymin}\n")
        f.write(f"set ymax {ymax}\n")
        f.write(f"set zmin 0.0 \n")
        f.write(f"set zmax 1\n")
        f.write("set grid [pw::Grid getAll]\n")
        f.write('set cells [pw::Grid getCellsWithinBoundingBox $grid "$xmin $ymin $zmin $xmax $ymax $zmax"]\n')
        f.write('set zoneName "actuatorDiskZone"\n')
        f.write("set zone [pw::CellZone create $zoneName]\n")
        f.write("$zone addCells $cells\n")
        #f.write("puts 'Created cell zone "$zoneName" with [llength $cells] cells.'\n")



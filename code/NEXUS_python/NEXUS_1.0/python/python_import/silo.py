

import numpy as np
import sys
import math
from scipy.weave import inline, converters
from miscellaneous import throwError, throwWarning

SILO_PATH = "/net/plato/data/users/cautun/Programs/stow/"


#~~~~~~~ Functions ~~~~~~~~

# writeHalos(fileName, pos=None,posUnit="Mpc", vel=None,velUnit="km/s", mass=None,massUnit="M0", rad=None,radUnit="Mpc", id=None, time=1.0, switchXYZ=True, VERBOSE=True)

# writeRectilinearMesh(fileName, meshName="mesh", posX, posY=None, posZ=None, posLabels=["X","Y","Z"], posUnits=["Mpc","Mpc","Mpc"], scalarData=None, scalarLabels=None, scalarDataUnits=None, vectorData=None, vectorLabels=None, vectorDataUnits=None, noVectorComponents=None, expressions=None, zoneCentered=True, VERBOSE=True)


SILO_VARTYPE = { "SCALAR":200, "VECTOR":201, "TENSOR":202, "SYMTENSOR":203, "ARRAY":204, "MATERIAL":205, "SPECIES":206, "LABEL":207 }
# inline(""" std::cout << DB_VARTYPE_SCALAR << "\\n"; """, [], headers=['<silo.h>'], include_dirs=[SILO_PATH+"include"], library_dirs=[SILO_PATH+"lib"], libraries=['m','silo'])


def writeHalos(fileName, pos=None,posUnit="Mpc", vel=None,velUnit="km/s", mass=None,massUnit="M0", rad=None,radUnit="Mpc", id=None, time=1.0, switchXYZ=True, VERBOSE=True):
    """ This function writes to a silo file the halo properties. """
    if VERBOSE:
        print "Writting halo data to the silo file '%s' ... " % fileName,
        sys.stdout.flush()
    
    # set some variables
    axisOrder = ["0","1","2"]
    if switchXYZ is True:
        axisOrder = ["2","1","0"]   # switches the X and Z axes - since the density reader does the same, we can compare easier between the two
    
    # use weave to write C code that will do the actual file writting
    siloWriteHaloSupportCode_header = """
    #line 1000 "siloWriteHaloSupportCode"
    #include <cstdio>
    #include <silo.h>
    #include <string>
    #include <vector>
    #include <new>
    """
    siloWriteHaloSupportCode_positions = """
    void haloPositions(DBfile *db, blitz::Array<float,2> pos, char const *posUnit, float time)  // this function writes the halo positions to a file
    {
        int noPoints = pos.extent(0);
        float *pointCoords[3];
        for (size_t i=0; i<3; ++i)
            pointCoords[i] = new float[noPoints];
        for (int i=0; i<noPoints; ++i)
        {
            pointCoords[0][i] = pos(i,%s);
            pointCoords[1][i] = pos(i,%s);
            pointCoords[2][i] = pos(i,%s);
        }
        
        // Write the point mesh
        DBoptlist *optlist = DBMakeOptlist(7);
        DBAddOption(optlist, DBOPT_XLABEL, (void*)"X");
        DBAddOption(optlist, DBOPT_YLABEL, (void*)"Y");
        DBAddOption(optlist, DBOPT_ZLABEL, (void*)"Z");
        DBAddOption(optlist, DBOPT_XUNITS, (void *)posUnit);
        DBAddOption(optlist, DBOPT_YUNITS, (void *)posUnit);
        DBAddOption(optlist, DBOPT_ZUNITS, (void *)posUnit);
        DBAddOption(optlist, DBOPT_TIME, (void *)(&time));
        DBPutPointmesh(db, "HaloMesh", 3, pointCoords, noPoints, DB_FLOAT, optlist);
        DBFreeOptlist(optlist);
        for (size_t i=0; i<3; ++i)
            delete[] pointCoords[i];
    }
    """ % (axisOrder[0], axisOrder[1], axisOrder[2])
    siloWriteHaloSupportCode_variable = """
    void haloVariable(DBfile *db,blitz::Array<float,1> var, char const *varUnit, float time, char *varName)  // this function writes a halo vector property to file
    {
        int noPoints = var.extent(0);
        float *varCoords[1];
        varCoords[0] = new float[noPoints];
        
        for (int i=0; i<noPoints; ++i)
            varCoords[0][i] = var(i);
        
        // Write the vector variable
        DBoptlist *optlist = DBMakeOptlist(2);
        DBAddOption(optlist, DBOPT_UNITS, (void *)varUnit);
        DBAddOption(optlist, DBOPT_TIME, (void *)(&time));
        DBPutPointvar(db, varName, "HaloMesh", 1, varCoords, noPoints, DB_FLOAT, optlist);
        DBFreeOptlist(optlist);
        delete[] varCoords[0];
    }
    void haloVariable(DBfile *db,blitz::Array<int,1> var, char const *varUnit, float time, char *varName)  // this function writes a halo vector property to file
    {
        int noPoints = var.extent(0);
        int *varCoords[1];
        varCoords[0] = new int[noPoints];
        
        for (int i=0; i<noPoints; ++i)
            varCoords[0][i] = var(i);
        
        // Write the vector variable
        DBoptlist *optlist = DBMakeOptlist(2);
        DBAddOption(optlist, DBOPT_UNITS, (void *)varUnit);
        DBAddOption(optlist, DBOPT_TIME, (void *)(&time));
        DBPutPointvar(db, varName, "HaloMesh", 1, varCoords, noPoints, DB_INT, optlist);
        DBFreeOptlist(optlist);
        delete[] varCoords[0];
    }
    void haloVariable(DBfile *db,blitz::Array<float,2> var, int const N, char const *varUnit, float time, char *varName)  // this function writes a halo vector property to file
    {
        int noPoints = var.extent(0);
        float *varCoords[N];
        for (int i=0; i<N; ++i)
            varCoords[i] = new float[noPoints];
        if (N==3)
            for (int i=0; i<noPoints; ++i)
            {
                varCoords[0][i] = var(i,%s);
                varCoords[1][i] = var(i,%s);
                varCoords[2][i] = var(i,%s);
            }
        else
            for (int i=0; i<noPoints; ++i)
                for (int j=0; j<N; ++j)
                    varCoords[j][i] = var(i,j);
        
        // Write the vector variable
        DBoptlist *optlist = DBMakeOptlist(2);
        DBAddOption(optlist, DBOPT_UNITS, (void *)varUnit);
        DBAddOption(optlist, DBOPT_TIME, (void *)(&time));
        DBPutPointvar(db, varName, "HaloMesh", N, varCoords, noPoints, DB_FLOAT, optlist);
        DBFreeOptlist(optlist);
        for (int i=0; i<N; ++i)
            delete[] varCoords[i];
    }
    """ % (axisOrder[0], axisOrder[1], axisOrder[2])
    siloWriteHaloSupportCode_expressions = """
    void haloExpressions(DBfile * db )
    {
        /* Write some expressions to the Silo file. */
        char *names[] = {"haloPos_X", "haloPos_Y", "haloPos_Z"};
        char *defs[] = {"coord(HaloMesh)[0]", "coord(HaloMesh)[1]", "coord(HaloMesh)[2]"};
        int types[] = {DB_VARTYPE_SCALAR, DB_VARTYPE_SCALAR, DB_VARTYPE_SCALAR};
        DBPutDefvars( db, "halo_position_variables", 3, names, types, defs, NULL );
    }
    """
    siloWriteHaloSupportCode = siloWriteHaloSupportCode_header + siloWriteHaloSupportCode_positions + siloWriteHaloSupportCode_variable + siloWriteHaloSupportCode_expressions
    
    siloWriteHaloCode = """
    #line 1000 "siloWriteHaloCode"
    DBfile *db;
    char const *file = fileName.c_str();
    
    db = DBCreate( file, DB_CLOBBER, DB_LOCAL, "Halo properties storage", DB_PDB);
    haloPositions( db, pos, posUnit.c_str(), float(time) );                           // write the halo positions
    if (velOn==1) haloVariable( db, vel, 3, velUnit.c_str(), float(time), (char *)"velocity" );  // write the halo velocities
    if (massOn==1) haloVariable( db, mass, massUnit.c_str(), float(time), (char *)"mass" );   // write the halo masses
    if (radOn==1) haloVariable( db, rad, radUnit.c_str(), float(time), (char *)"radius" );    // write the halo radius
    if (idOn==1) haloVariable( db, id, (char const *)"none", float(time), (char *)"id" );     // write the halo id
    haloExpressions(db);        // write some expressions for VisIt
    DBClose(db); 
    """
    
    # check which variables were supplied to the function
    if pos is None:
        throwError( "The function 'writeHalos' needs the halo positions to be able to write the halo properties in a file to be visualized with VisIt." )
    velOn, massOn, radOn, idOn = 0, 0, 0, 0
    if vel is not None: velOn = 1
    if mass is not None: massOn = 1
    if rad is not None: radOn = 1
    if id is not None:
        idOn = 1
    
    # write the results to file
    inline( siloWriteHaloCode, ['fileName', 'pos','posUnit', 'vel','velUnit','velOn', 'mass','massUnit','massOn', 'rad','radUnit','radOn', 'id','idOn', 'time'], type_converters=converters.blitz, support_code=siloWriteHaloSupportCode, include_dirs=[SILO_PATH+"include"], library_dirs=[SILO_PATH+"lib"], libraries=['m','silo'] )
    
    if VERBOSE: print "Done."


def writeRectilinearMesh(fileName, meshName="mesh", posX=None, posY=None, posZ=None, posLabels=["X","Y","Z"], posUnits=["Mpc","Mpc","Mpc"], scalarData=None, scalarLabels=None, scalarDataUnits=None, vectorData=None, vectorLabels=None, vectorDataUnits=None, noVectorComponents=None, expressions=None, zoneCentered=True, time=1.0, VERBOSE=True):
    """ Writes in a silo file a rectiliniar mesh. It takes the following argumnets (besides the obvious ones):
    meshName = name of the mesh - a string variable (no fancy symbols or space)
    posX, posY, posZ = give the coordinates of the nodes along each axis (numpy matrices) - set posY=None for outputing a 1D result and posZ=None for a 2D result
    posLabels = positions labels for the coordinates
    posUnits = units for the positions
    scalarData = the numpy matrix giving the scalar data to be written to file (can contain more than one variable - than it should have one additional dimensions on top of the spatial dimensions).
    scalarLabels = name of the scalar variables
    scalarDataUnits = if supplied, must have the same length as the number of scalar variables. A list of strings giving the unit for each scalar variable.
    vectorData = the numpy matrix giving the vector data to be written to file (can contain more than one variable - than it should have one additional dimensions on top of the spatial dimensions). By default, each vector has the same number of components as the spatial dimensions, unless specified otherwise in the 'noVectorComponents' variable.
    vectorLabels = name of the vector variables
    vectorDataUnits = see 'scalarDataUnits'.
    noVectorComponents = how many components are in a vector (if different from the number of spatial dimensions).
    expressions = a 2D numpy array giving the expressions for VisIt variables. each row should contain: expressionName expressionDefinition and expressionType (DB_VARTYPE_SCALAR, DB_VARTYPE_VECTOR, DB_VARTYPE_TENSOR, DB_VARTYPE_SYMTENSOR, DB_VARTYPE_ARRAY, DB_VARTYPE_MATERIAL, DB_VARTYPE_SPECIES, DB_VARTYPE_LABEL).
    zoneCentered = true (DEFAULT) if the variable is zone centered vs. node centered. Zone centered variables should have a (Nx-1)*(Ny-1)*(Nz-1) size versus the node centered variables which should have a Nx*Ny*Nz size.
     """
    
    # this is the C++ code that does the actual writing to file of the results
    writeRectilinearMeshSupportCode_header = """
    #line 1000 "writeRectilinearMeshSupportCode"
    #include <cstdio>
    #include <silo.h>
    #include <string>
    #include <vector>
    #include <new>
    """
    writeRectilinearMeshSupportCode_mesh = """
    void writeMesh(DBfile *db, char * mesh_name, int const N, blitz::Array<float,1> **coords,  char const **coords_label, char const **coords_unit, float time)  // this function writes the mesh to a file
    {
        int size[N];
        float *mesh[N];
        for (int i=0; i<N; ++i)
        {
            size[i] = coords[i]->extent(0);
            mesh[i] = new float[size[i]];
            for (int j=0; j<size[i]; ++j)
                mesh[i][j] = (*coords[i])(j);
        }
        
        // Write the mesh
        DBoptlist *optlist = DBMakeOptlist(2*N+1);
        DBAddOption(optlist, DBOPT_XLABEL, (void*)coords_label[0]);
        DBAddOption(optlist, DBOPT_XUNITS, (void *)coords_unit[0]);
        if ( N>1 )
        {
            DBAddOption(optlist, DBOPT_YLABEL, (void*)coords_label[1]);
            DBAddOption(optlist, DBOPT_YUNITS, (void *)coords_unit[1]);
        }
        if ( N>2 )
        {
            DBAddOption(optlist, DBOPT_ZLABEL, (void*)coords_label[2]);
            DBAddOption(optlist, DBOPT_ZUNITS, (void *)coords_unit[2]);
        }
        DBAddOption(optlist, DBOPT_TIME, (void *)(&time));
        DBPutQuadmesh(db, mesh_name, NULL, mesh, size, N, DB_FLOAT, DB_COLLINEAR, optlist);
        DBFreeOptlist(optlist);
        for (size_t i=0; i<N; ++i)
            delete[] mesh[i];
    }
    """
    writeRectilinearMeshSupportCode_variables = """
    void writeVariable(DBfile *db, char * mesh_name, int *dims, int const N, float **scalar, int const no_components, char *scalar_label, char const *scalar_unit, float time, int centering)  // this function writes the mesh to a file
    {
        // Write the variable onto the mesh
        char *varNames[no_components];
        for (int i=0; i<no_components; ++i)
        {
            varNames[i] = new char[100];
            std::sprintf(varNames[i], "%s_%i", scalar_label, i);
        }
        DBoptlist *optlist = DBMakeOptlist(2);
        DBAddOption(optlist, DBOPT_UNITS, (void *)scalar_unit);
        DBAddOption(optlist, DBOPT_TIME, (void *)(&time));
        DBPutQuadvar(db, scalar_label, mesh_name, no_components, varNames, scalar, dims, N, NULL, 0, DB_FLOAT, centering, optlist);
        DBFreeOptlist(optlist);
        for (int i=0; i<no_components; ++i) delete varNames[i];
    }
    /* Copies the 'component' of the last component in the input array. */
    template <typename T, int N>
    void copyData(blitz::Array<T,N> &data, int const component, float *output)
    {
        int const i = component;
        if (N==2)
            for (int i1=0; i1<data.extent(0); ++i1)
                output[i1] = data(i1,i);
        if (N==3)
        {
            int ny = data.extent(1);
            for (int i1=0; i1<data.extent(0); ++i1)
                for (int i2=0; i2<data.extent(1); ++i2)
                    output[i1*ny+i2] = data(i1,i2,i);
        }
        if (N==4)
        {
            int nz = data.extent(2), nyz = data.extent(1) * nz;
            for (int i1=0; i1<data.extent(0); ++i1)
                for (int i2=0; i2<data.extent(0); ++i2)
                    for (int i3=0; i3<data.extent(0); ++i3)
                        output[i1*nyz+i2*nz+i3] = data(i1,i2,i3,i);
        }
    }
    template <typename T, int N>
    void copyData(blitz::Array<T,N> &data, float *output)
    {
        if (N==1)
            for (int i1=0; i1<data.extent(0); ++i1)
                output[i1] = data(i1);
        if (N==2)
        {
            int ny = data.extent(1);
            for (int i1=0; i1<data.extent(0); ++i1)
                for (int i2=0; i2<data.extent(1); ++i2)
                    output[i1*ny+i2] = data(i1,i2);
        }
        if (N==3)
        {
            int nz = data.extent(2), nyz = data.extent(1) * nz;
            for (int i1=0; i1<data.extent(0); ++i1)
                for (int i2=0; i2<data.extent(0); ++i2)
                    for (int i3=0; i3<data.extent(0); ++i3)
                        output[i1*nyz+i2*nz+i3] = data(i1,i2,i3);
        }
    }
    void assignMemory(float **output, int const noComponents, int const size)
    {
        for (int i=0; i<noComponents; ++i)
            output[i] = new float[size];
    }
    void deleteMemory(float **output, int const noComponents)
    {
        for (int i=0; i<noComponents; ++i)
            delete[] output[i];
    }
    """
    writeRectilinearMeshSupportCode_expressions = """
    template <typename T>
    void writeExpressions(DBfile *db, T expressions, int const noExpressions)
    {
        /* Write some expressions to the Silo file. */
        char *names[noExpressions];
        char *defs[noExpressions];
        int types[noExpressions];
        for (int i=0; i<noExpressions; ++i)
        {
            PyObject *item = PyList_GetItem(expressions,i);
            names[i] = PyString_AsString(PyList_GetItem(item,0));
            defs[i] = PyString_AsString(PyList_GetItem(item,1));
            types[i] = atoi( PyString_AsString(PyList_GetItem(item,2)) );
        }
        DBPutDefvars( db, "data expressions", noExpressions, names, types, defs, NULL);
    }
    """
    writeRectilinearMeshSupportCode = writeRectilinearMeshSupportCode_header + writeRectilinearMeshSupportCode_mesh + writeRectilinearMeshSupportCode_variables + writeRectilinearMeshSupportCode_expressions
     
    writeRectilinearMeshCode = """
    #line 1000 "writeRectilinearMeshCode"
    // define some constants
    blitz::Array<float,1> *coords[3] = { &posX, &posY, &posZ };
    char const *coords_label[3] = { PyString_AsString(PyList_GetItem(posLabels,0)), PyString_AsString(PyList_GetItem(posLabels,1)), PyString_AsString(PyList_GetItem(posLabels,2)) };
    char const *coords_unit[3] = { PyString_AsString(PyList_GetItem(posUnits,0)), PyString_AsString(PyList_GetItem(posUnits,1)), PyString_AsString(PyList_GetItem(posUnits,2)) };
    
    int centering = DB_NODECENT, zonal = outputDataProperties(8);
    if (zonal==1) centering = DB_ZONECENT;
    
    int const noDims = outputDataProperties(0);
    int dims[noDims];
    int dataSize = 1;
    for (int i=0; i<noDims; ++i)
    {
        if (zonal==1) dims[i] = outputDataProperties(i+1) - 1;
        else dims[i] = outputDataProperties(i+1);
        dataSize *= dims[i];
    }
    int noScalarFields = outputDataProperties(4), noVecComp = outputDataProperties(5), noVectorFields = outputDataProperties(6), noExpressions = outputDataProperties(7);
    
    // open the output file
    DBfile *db;
    char const *file = fileName.c_str();
    
    db = DBCreate( file, DB_CLOBBER, DB_LOCAL, "Data mesh storage", DB_PDB);
    writeMesh( db, (char *)meshName.c_str(), noDims, coords, coords_label, coords_unit, float(time) );   // write the mesh coordinates
    
    // write the scalar variables
    for (int i=0; i<noScalarFields; ++i)
    {
        float *output[1];
        assignMemory( output, 1, dataSize );
        if (scalarData.dimensions()==noDims) copyData( scalarData, output[0] );
        else copyData( scalarData, i, output[0] );
        writeVariable( db, (char *)meshName.c_str(), dims, noDims, output, 1, PyString_AsString(PyList_GetItem(scalarLabels,i)), PyString_AsString(PyList_GetItem(scalarDataUnits,i)), float(time), centering );
        deleteMemory( output, 1 );
    }
    
    // write the vector variables
    for (int i=0; i<noVectorFields; ++i)
    {
        float *output[noVecComp];
        assignMemory( output, noVecComp, dataSize );
        for (int j=0; j<noVecComp; ++j)
            copyData( vectorData, i*noVecComp+j, output[j] );
        writeVariable( db, (char *)meshName.c_str(), dims, noDims, output, noVecComp, PyString_AsString(PyList_GetItem(vectorLabels,i)), PyString_AsString(PyList_GetItem(vectorDataUnits,i)), float(time), centering );
        deleteMemory( output, noVecComp );
    }
    
    // write some expressions for VisIt - bla bla
    if (noExpressions!=0)
        writeExpressions( db, expressions, noExpressions );
    DBClose(db);
    """
    
    
    
    if VERBOSE:
        print "Writing the data on a rectiliniar mesh to the silo file '%s': " % fileName
    
    # find the number of dimensions of the grid
    noDims, xSize, ySize, zSize = -1, 1, 1, 1
    if posX is None:
        throwError( "The function 'writeRectilinearMesh' needs at least the x-coordinates of the 1D mesh that will be written to the silo file." )
    elif posY is None:
        noDims = 1
        xSize = posX.size
        posY, posZ = np.zeros(2,np.float32), np.zeros(2,np.float32)
    elif posZ is None:
        noDims = 2
        xSize, ySize = posX.size, posY.size
        posZ = np.zeros(2,np.float32)
    else:
        noDims = 3
        xSize, ySize, zSize = posX.size, posY.size, posZ.size
    if VERBOSE: print "\t Found a %iD mesh of dimensions (%i x %i x %i) that needs to be written to file." % (noDims,xSize,ySize,zSize)
    # Check that everything is as expected
    if len(posX.shape)!=1: throwError( "The 'posX' array giving the grid coordinates of the x-axis needs to be a 1D numpy array." )
    if posY is not None and len(posY.shape)!=1: throwError( "The 'posY' array giving the grid coordinates of the y-axis needs to be a 1D numpy array." )
    if posZ is not None and len(posZ.shape)!=1: throwError( "The 'posZ' array giving the grid coordinates of the z-axis needs to be a 1D numpy array." )
    # make sure that the posLabels and posUnits have at least 3 elements
    if posX.dtype is not np.float32: posX = posX.astype( np.float32 )
    if posY.dtype is not np.float32: posY = posY.astype( np.float32 )
    if posZ.dtype is not np.float32: posZ = posZ.astype( np.float32 )
    coordinateNames = {0:"X",1:"Y",2:"Z"}
    for i in range(len(posLabels),3): posLabels.append( coordinateNames[i] )
    for i in range(len(posUnits),3): posUnits.append( "Mpc" )
    
    # check the scalar data
    expectedSize = xSize * ySize * zSize
    if zoneCentered:
        if noDims==1: expectedSize = xSize-1
        if noDims==2: expectedSize = (xSize-1) * (ySize-1)
        if noDims==3: expectedSize = (xSize-1) * (ySize-1) * (zSize-1)
    if expectedSize==0: throwError( "The expected data size is 0. This represents an error. Possible causes: one of the arrays giving the axis coordinates has a length of 1 and the data is zone centered." )
    noScalarFields = 0
    hasScalarUnits = False
    if scalarData is not None:
        if len(scalarData.shape)!=noDims and len(scalarData.shape)!=noDims+1:
            throwError( "The scalar data supplied to the function 'writeRectilinearMesh' in the 'scalarData' variable needs to be a %i dimensional array for a single scalar field or a %i dimensional array for multiple scalar field components." % (noDims,noDims+1) )
        elif len(scalarData.shape)==noDims:
            noScalarFields = 1
        else:
            noScalarFields = scalarData.shape[noDims]
        if scalarDataUnits is not None:
            if len(scalarDataUnits)!=noScalarFields:
                print "The 'scalarDataUnits' variable keeping track of the units for the scalar data must be the same as the number of scalar components. The output file will not contain units for the scalar fields since the 'scalarDataUnits' variable has a length of %i while there are %i scalar fields." % (len(scalarDataUnits),noScalarFields)
            else:
                hasScalarUnits = True
        scalarSize = scalarData.size / noScalarFields
        if scalarSize!=expectedSize:
            throwError( "The size of the scalar fields grid is %i while the expected size (using the input coordinate arrays) is %i. The program cannot continue due to the mismatch." % (scalarSize,expectedSize) )
        if hasScalarUnits and VERBOSE: print "\t Found %i scalar fields with units." % noScalarFields
        elif VERBOSE: print "\t Found %i scalar fields (NO units available)." % noScalarFields
        # make sure that scalarLabels and scalarDataUnits have entries for all scalar fields
        for i in range(len(scalarLabels),noScalarFields): scalarLabels.append( "scalar"+"%i" % i )
        for i in range(len(scalarDataUnits),noScalarFields): scalarDataUnits.append( "" )
    else: scalarData = np.zeros(2,np.float32)
    
    
    # check the vector data
    noVectorFields = 0
    hasVectorUnits = False
    if noVectorComponents is None: noVectorComponents = noDims
    if vectorData is not None:
        if len(vectorData.shape)!=noDims+1:
            throwError( "The vector data supplied to the function 'writeRectilinearMesh' in the 'vectorData' variable needs to be a %i dimensional array (if multiple vector fields are present, all their components must be written in the last dimension)." % noDims+1 )
        else:
            noVectorFields = vectorData.shape[noDims] / noVectorComponents
        if vectorDataUnits is not None:
            if len(vectorDataUnits)!=noVectorFields:
                print "The 'vectorDataUnits' variable keeping track of the units for the vector data must be the same as the number of vector components. The output file will not contain units for the vector fields since the 'vectorDataUnits' variable has a length of %i while there are %i vector fields." % (len(vectorDataUnits),noVectorFields)
            else:
                hasVectorUnits = True
        vectorSize = vectorData.size / (noVectorComponents* noVectorFields)
        if vectorSize!=expectedSize:
            throwError( "The size of the vector fields grid is %i while the expected size (using the input coordinate arrays) is %i. The program cannot continue due to the mismatch." % (vectorSize,expectedSize) )
        if hasVectorUnits and VERBOSE: print "\t Found %i vector fields of %i components with units." % (noVectorFields,noVectorComponents)
        elif VERBOSE: print "\t Found %i vector fields of %i components (NO units available)." % (noVectorFields,noVectorComponents)
        # make sure that vectorLabels and vectorDataUnits have entries for all vector fields
        for i in range(len(vectorLabels),noVectorFields): vectorLabels.append( "vector"+"%i" % i )
        for i in range(len(vectorDataUnits),noVectorFields): vectorDataUnits.append( "" )
    else: vectorData = np.zeros(2,np.float32)
    
    
    # check the expressions
    noExpressions = 0
    if expressions is not None:
        if len(expressions.shape)!=2:
            throwError( "The 'expressions' array supplied to the 'writeRectilinearMesh' function must be a 2D numpy array." )
        elif expressions.shape[1]!=3:
            throwError( "The 'expressions' array supplied to the 'writeRectilinearMesh' function must have 3 columns that give: expressionName, expressionDefinition and expressionType - in this order." )
        else:
            noExpressions = expressions.size / 3
            if VERBOSE: print "\t Found %i expressions." % noExpressions
        for i in range(noExpressions):
            if expressions[i,2] in SILO_VARTYPE.keys(): expressions[i,2] = "%i" % SILO_VARTYPE[expressions[i,2]]
            else: throwError( "Unregonized variable type '%s' in expression %i" % (expressions[i,2],i) )
    
    
    # write the results to file
    if VERBOSE: print "\t Writing the data to the silo file ...",
    zonal = 0
    if zoneCentered: zonal = 1
    outputDataProperties = np.array( [noDims,xSize,ySize,zSize, noScalarFields, noVectorComponents,noVectorFields, noExpressions, zonal], np.int32)
    inline( writeRectilinearMeshCode, ['fileName', 'meshName', 'outputDataProperties', 'posX', 'posY', 'posZ', 'posLabels', 'posUnits', 'scalarData', 'scalarLabels', 'scalarDataUnits', 'vectorData', 'vectorLabels', 'vectorDataUnits', 'expressions', 'time'], type_converters=converters.blitz, support_code=writeRectilinearMeshSupportCode, include_dirs=[SILO_PATH+"include"], library_dirs=[SILO_PATH+"lib"], libraries=['m','silo'] )
    
    
    if VERBOSE: print "Done."
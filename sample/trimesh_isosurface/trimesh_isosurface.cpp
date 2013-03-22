/****************************************************************************
* VCGLib o o *
* Visual and Computer Graphics Library o o *
* _ O _ *
* Copyright(C) 2004-2012 \/)\/ *
* Visual Computing Lab /\/| *
* ISTI - Italian National Research Council | *
* \ *
* All rights reserved. *
* *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 2 of the License, or *
* (at your option) any later version. *
* *
* This program is distributed in the hope that it will be useful, *
* but WITHOUT ANY WARRANTY; without even the implied warranty of *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt) *
* for more details. *
* *
****************************************************************************/
#include <limits>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <vcg/complex/complex.h>
#include <vcg/math/perlin_noise.h>
#include <vcg/complex/algorithms/create/marching_cubes.h>
//#include <vcg/complex/algorithms/create/extended_marching_cubes.h>
#include <vcg/complex/algorithms/create/mc_trivial_walker.h>
#define M_PI_2 1.57079632679489661923
//#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/smooth.h>

#include <GL/glut.h>
#include <wrap/io_trimesh/export_dae.h>
#include <wrap/io_trimesh/export_ply.h>

// update
#include <vcg/complex/algorithms/update/topology.h>

// local optimization
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

using namespace std;
using namespace vcg;

typedef float ScalarType;


class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>	::AsVertexType,
Use<MyFace>	::AsFaceType>{};

class MyVertex : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags , vertex::Normal3f, vertex::VFAdj, vertex::Mark, vertex::Color4b>{
      public:
        vcg::math::Quadric<double> &Qd() {return q;}
      private:
        math::Quadric<double> q;
};
class MyFace : public Face< MyUsedTypes, face::VertexRef, face::BitFlags, face::Normal3f, face::VFAdj , face::FFAdj> {};

class MyMesh	: public vcg::tri::TriMesh< std::vector< MyVertex>, std::vector< MyFace > > {};

typedef MyMesh::FaceIterator FaceIterator;
typedef MyMesh::FacePointer FacePointer;
typedef MyMesh::VertexPointer VertexPointer;
typedef MyMesh::VertexIterator VertexIterator;

namespace vcg
{
namespace tri
{

typedef BasicVertexPair<MyVertex> VertexPair;


typedef	SimpleTempData<MyMesh::VertContainer, math::Quadric<double> > QuadricTemp;


class QHelper
        {
        public:
      QHelper(){}
      static void Init(){}
      static math::Quadric<double> &Qd(MyVertex &v) {return TD()[v];}
      static math::Quadric<double> &Qd(MyVertex *v) {return TD()[*v];}
      static MyVertex::ScalarType W(MyVertex * /*v*/) {return 1.0;}
      static MyVertex::ScalarType W(MyVertex & /*v*/) {return 1.0;}
      static void Merge(MyVertex & /*v_dest*/, MyVertex const & /*v_del*/){}
      static QuadricTemp* &TDp() {static QuadricTemp *td; return td;}
      static QuadricTemp &TD() {return *TDp();}
        };


class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric< MyMesh, VertexPair , MyTriEdgeCollapse, QHelper > {
                        public:
            typedef vcg::tri::TriEdgeCollapseQuadric< MyMesh, VertexPair, MyTriEdgeCollapse, QHelper> TECQ;
            inline MyTriEdgeCollapse( const VertexPair &p, int i, BaseParameterClass *pp) :TECQ(p,i,pp){}
};
}
}
typedef SimpleVolume<SimpleVoxel> MyVolume;

class DAEConverter{ //reference vcg::tri::exportPLY
public:
    void saveModel(const char* file_name , MyMesh &_mesh  ,  void(*func)(int)){

        //_mesh = MyMesh(mesh);
        if( _mesh.FN() && _mesh.VN() ){
            SimpleTempData<MyMesh::VertContainer,int> indices(_mesh.vert);

            std::fstream fout;
            fout.open( file_name , std::fstream::out | std::fstream::trunc );
                fout << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" <<
                        "<COLLADA xmlns=\"http://www.collada.org/2005/11/COLLADASchema\" version=\"1.4.1\">\n" <<
                        "<library_visual_scenes>\n" <<
                        "<visual_scene id=\"ID1\">\n" <<
                        "<node name=\"teeth\">\n" <<
                        "<instance_geometry url=\"#ID2\">\n" <<
                        "</instance_geometry>\n" <<
                        "</node>\n" <<
                        "</visual_scene>\n" <<
                        "</library_visual_scenes>\n" <<
                        "<library_geometries>\n<geometry id=\"ID2\">\n<mesh>" << std::endl;
                func(10);
                this->writeSourceBlock(fout,_mesh ,indices);
                func(30);
                this->writeVertexBlock(fout,_mesh);
                func(60);
                this->writeTriangleBlock(fout,_mesh ,indices);
                func(90);
                fout << "</mesh>\n</geometry>\n" <<
                        "</library_geometries>\n" <<
                        "<scene>\n" <<
                        "<instance_visual_scene url=\"#ID1\" />" <<
                        "</scene>\n" <<
                        "</COLLADA>" << std::endl;
                fout.close();
                func(100);

        }
    }
    void writeSourceBlock( std::fstream & fout , MyMesh &_mesh, SimpleTempData<MyMesh::VertContainer,int> &indices){

        //vertexs
        fout << "<source id=\"ID5\">" <<
                "<float_array id=\"ID7" << "\" count=\"" << 3*_mesh.VN() << "\">";
        //write points


        VertexPointer  vp;
        VertexIterator vi;
        int j;
        for(j=0,vi=_mesh.vert.begin();vi!=_mesh.vert.end();++vi,j++){
            vp=&(*vi);
            indices[vi] = j;
            fout << setprecision(5) << vp->P()[0] << " " << setprecision(5) <<  vp->P()[1] << " " << setprecision(5) << vp->P()[2] << " ";
        }

        fout <<"</float_array>\n"
             << "<technique_common>\n" <<
                "<accessor count=\"" << _mesh.VN() << "\" source=\"#ID7" << "\" stride=\"" << 3 << "\">" << std::endl;
        fout << "<param name=\"X\" type=\"float\" />\n" <<
                "<param name=\"Y\" type=\"float\" />\n" <<
                "<param name=\"Z\" type=\"float\" />\n" <<
        "</accessor>\n</technique_common>\n</source>" << std::endl;
        //normals
        fout << "<source id=\"ID6\">" <<
                "<float_array id=\"ID8" << "\" count=\"" << 3*_mesh.FN() << "\">";
        //write normals
        FacePointer fp;

        FaceIterator fi;
        for(j=0,fi=_mesh.face.begin();fi!=_mesh.face.end();++fi){
            fp=&(*fi);
            fout << setprecision(5) << (double)fp->N()[0] << " " << setprecision(5) << (double)fp->N()[1] << " " << setprecision(5) << (double)fp->N()[2] << " ";
        }

        fout <<"</float_array>\n"
             << "<technique_common>\n" <<
                "<accessor count=\"" << _mesh.FN() << "\" source=\"#ID8" << "\" stride=\"" << 3 << "\">" << std::endl;
        fout << "<param name=\"X\" type=\"float\" />\n" <<
                "<param name=\"Y\" type=\"float\" />\n" <<
                "<param name=\"Z\" type=\"float\" />\n" <<
        "</accessor>\n</technique_common>\n</source>" << std::endl;
    }
    void writeVertexBlock( std::fstream &fout ,MyMesh &_mesh){
        fout << "<vertices id=\"ID9\">\n" <<
                "<input semantic=\"POSITION\" source=\"#ID5\" />\n" <<
                "<input semantic=\"NORMAL\" source=\"#ID6\" />\n"   <<
                "</vertices>" << std::endl;
    }
    void writeTriangleBlock( std::fstream &fout,MyMesh &_mesh , SimpleTempData<MyMesh::VertContainer,int> &indices ){
        FacePointer fp;
        int vv[3];
        FaceIterator fi;
        fout << "<triangles count=\"" << _mesh.FN() << "\">\n" <<
            "<input offset=\"0\" semantic=\"VERTEX\" source=\"#ID9\" />\n<p>";
        int j;
        for(j=0,fi=_mesh.face.begin();fi!=_mesh.face.end();++fi,j++){
            fp = &(*fi);
            vv[0]=indices[fp->cV(0)];
            vv[1]=indices[fp->cV(1)];
            vv[2]=indices[fp->cV(2)];
            fout << vv[0] << " " << vv[1] << " " << vv[2] << " ";
        }
        fout << "</p>\n</triangles>" << std::endl;

    }
};

void update(int v){
    printf("complete; %d\n" , v);
}

int main(int argc , char **argv)
{
MyVolume	volume;

  typedef vcg::tri::TrivialWalker<MyMesh,MyVolume>	MyWalker;
typedef vcg::tri::MarchingCubes<MyMesh, MyWalker>	MyMarchingCubes;
MyWalker walker;

    int w,h,d;
    float iso_value;
    const char* path;
    if( argc > 1){
        w = atoi( argv[1] );
        h = atoi( argv[2] );
        d = atoi( argv[3] );
        iso_value = atof( argv[4] );
        path = argv[5];
    }else{
        printf("unexpected parameters...");
        exit(1);
    }
  // Simple initialization of the volume with some cool perlin noise
    printf( "%d %d %d %.2f path: %s\n" , w, h, d, iso_value, path);
    volume.Init(Point3i(w,h,d));

    FILE* pfin = fopen( path , "rb");
    if( pfin == NULL )exit(1);

    float v;
    for(int k=0;k<d;k++)
        for(int j=0;j<h;j++)
            for(int i=0;i<w;i++){
                fread( &v , sizeof(float) , 1 , pfin );
            volume.Val(i,j,k)= 1 - v;
        }
    fclose( pfin );

// MARCHING CUBES
MyMesh	mc_mesh;
printf("[MARCHING CUBES] Building mesh...");
MyMarchingCubes	mc(mc_mesh, walker);
    walker.BuildMesh<MyMarchingCubes>(mc_mesh, volume, mc, iso_value);
    printf("generate face num %d, vertex num %d" , mc_mesh.vn , mc_mesh.fn );
    tri::UpdateTopology<MyMesh>::VertexFace(mc_mesh);
    tri::UpdateFlags<MyMesh>::FaceBorderFromVF(mc_mesh);
    math::Quadric<double> QZero;
    QZero.SetZero();
    tri::QuadricTemp TD(mc_mesh.vert,QZero);
    tri::QHelper::TDp()=&TD;

    TriEdgeCollapseQuadricParameter qparams;
    qparams.QualityThr =.3;
    //float TargetError=std::numeric_limits<float>::max();
    int FinalSize=mc_mesh.fn*0.2;
    printf("target face num %d\n" , FinalSize);
    qparams.PreserveTopology = true;
    //qparams.QuadricEpsilon=TargetError;
    qparams.NormalCheck=true;

    vcg::tri::UpdateBounding<MyMesh>::Box(mc_mesh);
    // decimator initialization
    vcg::LocalOptimization<MyMesh> DeciSession(mc_mesh,&qparams);

    int t1=clock();
    DeciSession.Init<MyTriEdgeCollapse>();
    int t2=clock();
    printf("Initial Heap Size %i\n",int(DeciSession.h.size()));

    DeciSession.SetTargetSimplices(FinalSize);
    DeciSession.SetTimeBudget(0.5f);

    //if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);

    printf("start mesh decimator....\n");
    while(DeciSession.DoOptimization() && mc_mesh.fn>FinalSize )
      printf("Current Mesh size %7i heap sz %9i err \r",mc_mesh.fn, int(DeciSession.h.size()));

    int t3=clock();
    printf("mesh %d %d Error %g \n",mc_mesh.vn,mc_mesh.fn,DeciSession.currMetric);
    printf("\nCompleted in (%i+%i) msec\n",t2-t1,t3-t2);
    DeciSession.Finalize<tri::MyTriEdgeCollapse >();
     tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(mc_mesh);

    printf("start taubin smooth....\n");
    //smooth
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(mc_mesh);
          int stepSmoothNum = 15;
          float lambda=0.5;
          float mu=-0.53;

    size_t cnt=tri::UpdateSelection<MyMesh>::VertexFromFaceStrict(mc_mesh);
    tri::Smooth<MyMesh>::VertexCoordTaubin(mc_mesh,stepSmoothNum,lambda,mu,cnt>0);
    //Log( "Smoothed %d vertices", cnt>0 ? cnt : mc_mesh.vn);
    //tri::UpdateTopology<MyMesh>::VertexFace(mc_mesh);
    tri::UpdateNormal<MyMesh>::PerFaceNormalized(mc_mesh);

    /* Explicitly delete unused vectors and faces,which is not impicitly cleared by exporterDAE */
    vcg::tri::Allocator<MyMesh>::CompactFaceVector(mc_mesh);
    vcg::tri::Allocator<MyMesh>::CompactVertexVector(mc_mesh);

    int mask = tri::io::Mask::IOM_FACENORMAL;
    /*
    vcg::tri::io::ExporterPLY<MyMesh>::Save( mc_mesh, "marching_cubes.dae" );

    DAEConverter dae_saver;

    dae_saver.saveModel( "my.dae" , mc_mesh , update );
    */
    vcg::tri::io::ExporterDAE<MyMesh>::Save( mc_mesh, "marching_cubes.dae" , mask );

printf("OK!\n");
}


/*!
   @file   SceneData.h
   @author Nong Li (Fall 2008)
   @author Remco Chang (Dec 2013)

   @brief  Header file containing scene data structures.
*/

#ifndef SCENEDATA_H
#define SCENEDATA_H

/* Includes */
#include <vector>
#include <string>
#include <glm/glm.hpp>
#include "../Shape.h"

using namespace std;

//! Enumeration for light types.
enum LightType {
   LIGHT_POINT, LIGHT_DIRECTIONAL, LIGHT_SPOT, LIGHT_AREA
};

//! Enumeration for types of primitives that can be stored in a scene file.
//enum PrimitiveType {
//	SHAPE_CUBE = 0,
//	SHAPE_CYLINDER = 1,
//	SHAPE_CONE = 2,
//	SHAPE_SPHERE = 3,
//	SHAPE_SPECIAL1 = 4,
//	SHAPE_SPECIAL2 = 5,
//	SHAPE_SPECIAL3 = 6,
//	SHAPE_MESH = 7
//};


//! Enumeration for types of transformations that can be applied to objects, lights, and cameras.
enum TransformationType {
   TRANSFORMATION_TRANSLATE, TRANSFORMATION_SCALE, 
   TRANSFORMATION_ROTATE, TRANSFORMATION_MATRIX
};

//! Struct to store a RGBA color in floats [0,1]
class SceneColor 
{
public:
    union {
        struct {
           float r;
           float g;
           float b;
           float a;
        };
        float channels[4]; // points to the same four floats above...
    };

   // @TODO: [OPTIONAL] You can implement some operators here for color arithmetic.

};

//! Scene global color coefficients
class SceneGlobalData 
{
public:
   float ka;  //! global ambient coefficient
   float kd;  //! global diffuse coefficient
   float ks;  //! global specular coefficient
   float kt;  //! global transparency coefficient
};

//! Data for a single light
class SceneLightData 
{
public:
   int id;
   LightType type;

   SceneColor color;
   glm::vec3 function;

   glm::vec3  pos;        //! Not applicable to directional lights
   glm::vec3 dir;         //! Not applicable to point lights

   float radius;        //! Only applicable to spot lights
   float penumbra;      //! Only applicable to spot lights
   float angle;         //! Only applicable to spot lights

   float width, height; //! Only applicable to area lights
};

//! Data for scene camera
class SceneCameraData
{
public:
   glm::vec3 pos;
   glm::vec3 lookAt;
   glm::vec3 look;
   glm::vec3 up;

   bool isDir;  //true if using look vector, false if using lookAt point

   float heightAngle;
   float aspectRatio;

   float aperture;      //! Only applicable for depth of field
   float focalLength;   //! Only applicable for depth of field
};

//! Data for file maps (ie: texture maps)
class SceneFileMap
{
public:
   bool isUsed;
   string filename;
   float repeatU;
   float repeatV;
};

//! Data for scene materials
class SceneMaterial 
{
public:
   SceneColor cDiffuse;
   SceneColor cAmbient;
   SceneColor cReflective;
   SceneColor cSpecular;
   SceneColor cTransparent;
   SceneColor cEmissive;

   SceneFileMap* textureMap;
   float blend;

   SceneFileMap* bumpMap;

   float shininess;

   float ior;           //! index of refaction
};

//! Data for a single primitive.
class ScenePrimitive
{
public:
   OBJ_TYPE type;
   string meshfile;     //! Only applicable to meshes
   SceneMaterial material;
};

/*!

@struct CS123SceneTransformation
@brief Data for transforming a scene object.

  Aside from the TransformationType, the remaining of the data in the
  struct is mutually exclusive

*/
class SceneTransformation
{
public:
   TransformationType type;

   /*! Translate type */
   glm::vec3 translate;

   /*! Scale type */
   glm::vec3 scale;

   /*! Rotate type */
   glm::vec3 rotate;
   float angle;        //! the angle of rotation, in **radians**

   /*! Matrix type */
   glm::mat4 matrix;
};

//! Structure for non-primitive scene objects
class SceneNode
{
public:

   /*! Transformation at this node */
   std::vector<SceneTransformation*> transformations;

   /*! Primitives at this node */
   std::vector<ScenePrimitive*> primitives;

   /*! Children of this node */
   std::vector<SceneNode*> children;
};

#endif


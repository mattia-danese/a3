/*!
   @file   SceneParser.h
   @author Eric Tamura (October 2006)
   @author Nong Li (August 2007)
   @author Remco Chang (December 2013)
*/

#ifndef SCENEPARSER_H
#define SCENEPARSER_H

#include "SceneData.h"

#include <vector>
#include <map>

class TiXmlElement;


//! This class parses the scene graph specified by the CS123 Xml file format.
class SceneParser
{
   public:
      //! Create a parser, passing it the scene file.
      SceneParser(const std::string& filename);
      //! Clean up all data for the scene
      ~SceneParser();

      //! Parse the scene.  Returns false if scene is invalid.
      bool parse();

      //! Returns global scene data
      bool getGlobalData(SceneGlobalData& data);

      //! Returns camera data
      bool getCameraData(SceneCameraData& data);

      //! Returns the root scenegraph node
      SceneNode* getRootNode();

      //! Returns number of lights in the scene
      int getNumLights();

      //! Returns the ith light data
      bool getLightData(int i, SceneLightData& data);


   private:
      //the filename should be contained within this parser implementation
      //if you want to parse a new file, instantiate a different parser
      bool parseGlobalData(const TiXmlElement* scenefile);
      bool parseCameraData(const TiXmlElement* scenefile);
      bool parseLightData(const TiXmlElement* scenefile);
      bool parseObjectData(const TiXmlElement* scenefile);
      bool parseTransBlock(const TiXmlElement* transblock, SceneNode* node);
      bool parsePrimitive(const TiXmlElement* prim, SceneNode* node);
      
      std::string file_name;
      mutable std::map<std::string, SceneNode*> m_objects;
      SceneCameraData m_cameraData;
      std::vector<SceneLightData*> m_lights;
      SceneGlobalData m_globalData;
      std::vector<SceneNode*> m_nodes;
};

#endif


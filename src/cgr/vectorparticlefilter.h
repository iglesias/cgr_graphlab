//========================================================================
//  This software is free: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License Version 3,
//  as published by the Free Software Foundation.
//
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  Version 3 in the file COPYING that came with this distribution.
//  If not, see <http://www.gnu.org/licenses/>.
//========================================================================
/*!
\file    vectorparticlefilter.h
\brief   C++ Interface: Particle, VectorLocalization
\author  Joydeep Biswas, (C) 2010-2012
*/
//========================================================================

#ifndef VECTORPARTICLEFILTER_H
#define VECTORPARTICLEFILTER_H

#include "stdio.h"
#include <graphlab.hpp>
#include <vector>
#include <map>
#include "vector_map.h"
#include <eigen3/Eigen/Dense>
#include "geometry.h"
#include "util.h"
#include "terminal_utils.h"

static const bool EnableProfiling = false;

//Forward declaration
class Particle2D;

typedef graphlab::distributed_graph<Particle2D, graphlab::empty> graph_type;

using namespace std;
using namespace Eigen;
/**
A Particle for 2D localization
**/
class Particle2D {
public:
  vector2f loc, lastLoc;
  float angle;
  float weight;
public:
  Particle2D() {weight = angle = 0.0; loc.zero();}
  Particle2D(float _x, float _y, float _theta, float _w) { loc.set(_x,_y); lastLoc.set(-DBL_MAX,-DBL_MAX); angle = _theta; weight = _w;}
  bool operator<(const Particle2D &other) {return weight<other.weight;}
  bool operator>(const Particle2D &other) {return weight>other.weight;}

  Particle2D &operator=(const Particle2D &other) {
    angle = other.angle;
    weight = other.weight;
    loc.x = other.loc.x;
    loc.y = other.loc.y;
    lastLoc.x = other.lastLoc.x;
    lastLoc.y = other.lastLoc.y;
    return *this;
  } 

  void save(graphlab::oarchive &oarc) const {
    oarc << angle << weight << loc.x << loc.y << lastLoc.x << lastLoc.y;
  }

  void load(graphlab::iarchive &iarc) {
    iarc >> angle >> weight >> loc.x >> loc.y >> lastLoc.x >> lastLoc.y;
  }
};

/**
Particle filter for vector localization
**/

class VectorLocalization2D{

public:
  typedef struct {
    double tStamp;

    //Parameters
    /// Alpha1 = Error in angle per delta in angle
    float Alpha1;
    /// Alpha2 = Error in angle per delta in translation
    float Alpha2;
    /// Alpha3 = Error in translation per delta in translation
    float Alpha3;

    float kernelSize;
  } MotionModelParams;

  class LidarParams{
    public:
    float* laserScan;
    int numRays;
    float minAngle;
    float maxAngle;
    float angleResolution;
    float minRange;
    float maxRange;
    int minPoints;
    Vector2f laserToBaseTrans;
    Matrix2f laserToBaseRot;
    vector<Vector2f> scanHeadings;
    
    int numSteps;
    float etaAngle;
    float etaLoc;
    float maxLocGradient;
    float maxAngleGradient;
    float attractorRange;
    float minCosAngleError;
    float correspondenceMargin;
    float minRefineFraction;
    
    float logObstacleProb; //Probability of an obstacle
    float logShortHitProb;
    float logOutOfRangeProb;
    float lidarStdDev;
    float correlationFactor;
    
    float kernelSize;
    
    void initialize();
  };
  
  class PointCloudParams{
    public: 
    double tStamp;
    float minRange;
    float maxRange;
    float fieldOfView;
    vector<vector2f> pointCloud;
    vector<vector2f> pointNormals;
    
    //Tunable parameters, MUST be set!
    float logObstacleProb; //Probability of an obstacle
    float logShortHitProb;
    float logOutOfRangeProb;
    float attractorRange;
    float stdDev;
    float corelationFactor;
    
    int numSteps;
    int minPoints;
    float etaAngle;
    float etaLoc;
    float maxLocGradient;
    float maxAngleGradient;
    float minCosAngleError;
    float correspondenceMargin;
    float minRefineFraction;
    
    float kernelSize;
  };
  
  
  typedef struct {
    double lastRunTime;
    double runTime;
    int numObservedPoints;
    int numCorrespondences;
    float stage0Weights;
    float stageRWeights;
    float meanSqError;
  } EvalValues;
  
  enum Resample{
    NaiveResampling,
    LowVarianceResampling,
    SensorResettingResampling,
  };

  struct ParticleInitializer {
    ParticleInitializer();
    ParticleInitializer(vector2f _loc, float _angle, float _locationUncertainty, float _angleUncertainty);

    vector2f loc;
    float angle;
    float locationUncertainty;
    float angleUncertainty;
  };

  struct Motion {
    Motion();
    Motion(float _dx, float _dy, float _dtheta, const MotionModelParams &_motionParams);

    double dx;
    double dy;
    double dtheta;
    const MotionModelParams* motionParams;
  };

  struct Refinement {
    Refinement();
    Refinement(const std::vector<vector2f>& _pointCloud, const std::vector<vector2f>& _pointNormals, const PointCloudParams& _pointCloudParams);

    const std::vector<vector2f>* pointCloud;
    const std::vector<vector2f>* pointNormals;
    const PointCloudParams* pointCloudParams;
  };

  struct Update {
    Update();
    Update(const MotionModelParams& _motionParams, const std::vector<vector2f>& _pointCloud, const std::vector<vector2f>& _pointNormals, const PointCloudParams& _pointCloudParams);

    const MotionModelParams* motionParams;
    const std::vector<vector2f>* pointCloud;
    const std::vector<vector2f>* pointNormals;
    const PointCloudParams* pointCloudParams;
  };

protected:
  //Distributed graph
  graph_type* graph;

  //Current state
  VectorMap* currentMap;
  vector2f currentLocation;
  float currentAngle;
  vector2f currentLocStdDev;
  float currentAngleStdDev;
  int numParticles;
  vector2f lastDistanceMoved;
  float lastAngleTurned;

  string mapsFolder;
  vector<VectorMap> maps;

  //These are class-wide only so that they can be accessed for debugging purposes
  vector<Vector2f> gradients;
  vector<Vector2f> points;

  //Statistics of performance
  int numUnrefinedParticlesSampled;
  int numRefinedParticlesSampled;
  float refinedImportanceWeights;
  float unrefinedImportanceWeights;
  double refineTime;
  double updateTime;
  EvalValues pointCloudEval;
  EvalValues laserEval;
    
public:
  VectorLocalization2D(int _numParticles, graph_type& _graph, const char* _mapsFolder);
  
  /// Sets Particle Filter LIDAR parameters
  void setParams(MotionModelParams _predictParams, LidarParams _lidarUpdateParams);
  /// Loads All the floor maps listed in atlas.txt
  void loadAtlas();
  /// Initialise arrays, and sets initial location to
  void initialize(const char* mapName, vector2f loc, float angle, float locationUncertainty = 0.0, float angleUncertainty = 0.0);
  /// Predict step of the particle filter. Samples from the motion model
  void predict(float dx, float dy, float dtheta, const VectorLocalization2D::MotionModelParams& motionParams);
  /// Refine proposal distribution based on LIDAR observations
  void refineLidar(const VectorLocalization2D::LidarParams& lidarParams);
  /// Refine proposal distribution based on Point Cloud observations
  void refinePointCloud(const vector<vector2f>& pointCloud, const vector< vector2f >& pointNormals, const VectorLocalization2D::PointCloudParams& pointCloudParams);
  /// Update distribution based on LIDAR observations
  void updateLidar(const VectorLocalization2D::LidarParams& lidarParams, const VectorLocalization2D::MotionModelParams& motionParams);
  /// Update distribution based on Point Cloud observations
  void updatePointCloud(const vector< vector2f >& pointCloud, vector< vector2f >& pointNormals, const VectorLocalization2D::MotionModelParams& motionParams, const VectorLocalization2D::PointCloudParams& pointCloudParams);
  /// Resample distribution
  void resample(Resample type = LowVarianceResampling);
  
  /// Refine a single location hypothesis based on a LIDAR observation
  void refineLocationLidar(vector2f& loc, float& angle, float& initialWeight, float& finalWeight, const VectorLocalization2D::LidarParams& lidarParams, const std::vector< Vector2f >& laserPoints);
  /// Refine a single location hypothesis based on a Point Cloud observation
  void refineLocationPointCloud(int particleidx, vector2f& loc, float& angle, float& initialWeight, float& finalWeight, const vector< vector2f >& pointCloud, const vector< vector2f >& pointNormals, const VectorLocalization2D::PointCloudParams& pointCloudParams) const;
  
  /// Attractor function used for refining location hypotheses 
  inline Vector2f attractorFunction(line2f l, Vector2f p, float attractorRange, float margin = 0) const;
  /// Observation function for a single ray
  inline Vector2f observationFunction(line2f l, Vector2f p) const;
  /// Gradient based on pointCloud observation
  void getPointCloudGradient(int particleIdx, vector2f loc, float angle, vector2f& locGrad, float& angleGrad, const std::vector< vector2f >& pointCloud, const std::vector< vector2f >& pointNormals, float& logWeight, const VectorLocalization2D::PointCloudParams& pointCloudParams, const std::vector< int >& lineCorrespondences, const std::vector< line2f >& lines) const;
  /// Gradient based on LIDAR observation
  void getLidarGradient(vector2f loc, float angle, vector2f& locGrad, float& angleGrad, float& logWeight, VectorLocalization2D::LidarParams lidarParams, const vector< Vector2f >& laserPoints, const vector<int> & lineCorrespondences, const vector<line2f> &lines);
  /// Observation likelihood based on LIDAR obhservation
  float observationWeightLidar(vector2f loc, float angle, const VectorLocalization2D::LidarParams& lidarParams, const std::vector< Vector2f >& laserPoints);
  /// Observation likelihood based on point cloud obhservation
  float observationWeightPointCloud(vector2f loc, float angle, const vector< vector2f >& pointCloud, const vector< vector2f >& pointNormals, const PointCloudParams& pointCloudParams) const;
  /// Probability of specified pose corresponding to motion model
  float motionModelWeight(vector2f loc, float angle, const VectorLocalization2D::MotionModelParams& motionParams) const;
  /// Set pose with specified uncertainty
  void setLocation(vector2f loc, float angle, const char* map, float locationUncertainty, float angleUncertainty);
  /// Set pose and map with specified uncertainty
  void setLocation(vector2f loc, float angle, float locationUncertainty, float angleUncertainty);
  /// Switch to a different map
  void setMap(const char * map);
  /// Resample particles using low variance resampling
  void lowVarianceResample();
  /// Resample particles using naive resampling
  void naiveResample();
  /// Compute the maximum likelihood location based on particle spread
  void computeLocation(vector2f &loc, float &angle);
  /// Returns the current map name
  const char* getCurrentMapName(){return currentMap->mapName.c_str();}
  /// Write to file run statistics about particle distribution
  void saveRunLog(FILE* f);
  /// Write to file riun-time profiling information
  void saveProfilingStats(FILE* f);
  /// Compile lists of drawing primitives that can be visualized for debugging purposes
  void drawDisplay(vector<float> &lines_p1x, vector<float> &lines_p1y, vector<float> &lines_p2x, vector<float> &lines_p2y, vector<uint32_t> &lines_color,
                   vector<float> &points_x, vector<float> &points_y, vector<uint32_t> &points_color, 
                   vector<float> &circles_x, vector<float> &circles_y, vector<uint32_t> &circles_color, float scale=1.0);
  /// Return evaluation values
  void getEvalValues(EvalValues &_laserEval, EvalValues &_pointCloudEval);
  /// Return angle and location uncertainties
  void getUncertainty(float &_angleUnc, float &_locUnc);
  /// Removes duplicate points with the same observation angle and range
  void reducePointCloud(const vector< vector2f >& pointCloud, const vector< vector2f >& pointNormals, vector< vector2f >& reducedPointCloud, vector< vector2f >& reducedPointNormals);
  /// Returns current particles
  void getParticles(std::vector<Particle2D> &_particles);
  /// Return the number of particles
  int getNumParticles() const { return numParticles; }
};

/// Initialize particle using location and uncertainties
void initializeParticle(graph_type::vertex_type& v);

/// Predict particle motion by sampling from the motion model
void predictParticle(graph_type::vertex_type& v);

/// Refine particle using point cloud observation
void refinePointCloudParticle(graph_type::vertex_type& v);

/// Compute particle sampling densities
void samplingDensityPointCloudParticle(graph_type::vertex_type& v);

/// Compute particle weights
void updatePointCloudParticle(graph_type::vertex_type& v);

/// MapReduce to compute the pose applying using maximum likelihood over the particle set
struct PoseReducer : public graphlab::IS_POD_TYPE {
  double x;
  double y;
  double heading_x;
  double heading_y;

  static PoseReducer getPose(const graph_type::vertex_type& v);

  PoseReducer& operator+=(const PoseReducer& other);
}; // struct PoseReducer

/// MapReduce to compute the total sampling density that will be used for weight normalization
struct SamplingDensityReducer : public graphlab::IS_POD_TYPE {
  float samplingDensity;

  static SamplingDensityReducer computeTotalDensity(const graph_type::vertex_type& v);

  SamplingDensityReducer& operator+=(const SamplingDensityReducer& other);
};

#endif //VECTORPARTICLEFILTER_H


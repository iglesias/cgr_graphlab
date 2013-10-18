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
\brief   C++ Implementation: Particle, VectorLocalization
\author  Joydeep Biswas, (C) 2010-2012
*/
//========================================================================

#include "vectorparticlefilter.h"

//Parameters required to initialize particles
const VectorLocalization2D::ParticleInitializer* PARTICLE_INITIALIZER = NULL;
//Motion currently processed
const VectorLocalization2D::Motion* MOTION = NULL;
//Current refinement parameters
const VectorLocalization2D::Refinement* REFINEMENT = NULL;
//Curent update parameters
const VectorLocalization2D::Update* UPDATE = NULL;
//Total sampling density used for normalization
float TOTAL_WEIGHT;

//This guy is created in localization_main
extern VectorLocalization2D* localization;

inline float eigenCross(const Vector2f &v1, const Vector2f &v2)
{
  return v1.dot(Vector2f(v2.y(),-v2.x()));
}

void VectorLocalization2D::LidarParams::initialize()
{
  float a0 = -0.5*angleResolution*float(numRays);
  scanHeadings.resize(numRays);
  int i=0;
  for(float a=minAngle; a<maxAngle; a+=angleResolution){
    scanHeadings[i] = Vector2f(cos(a),sin(a));
    i++;
  }
}

VectorLocalization2D::VectorLocalization2D(int _numParticles, graph_type& _graph, const char* _mapsFolder) {
  assert(_numParticles > 0);

  mapsFolder = string(_mapsFolder);
  loadAtlas();
  numParticles = _numParticles;

  //create distributed graph
  for (int id = 0; id < numParticles-1; ++id) {
    _graph.add_edge(id, id+1);
    _graph.add_edge(id+1, id);
  }

  //commit the distributed graph, denoting that it is no longer to be modified
  _graph.finalize();

  graph = &_graph;
  engine = new engine_type(_graph.dc(), _graph, "sync");
}

VectorLocalization2D::~VectorLocalization2D() {
  delete engine;
}

void VectorLocalization2D::loadAtlas()
{
  static const bool debug = true;
  string atlasFile = mapsFolder + "/atlas.txt";
  FILE* fid = fopen(atlasFile.c_str(),"r");
  if(fid==NULL){
    TerminalWarning("Unable to load Atlas!");
    return;
  }
  char mapName[4096];
  int mapNum;
  if(debug) printf("Loading Atlas...\n");
  maps.clear();
  while(fscanf(fid,"%d %s\n",&mapNum,mapName)==2){
    if(debug) printf("Loading map %s\n",mapName);
    maps.push_back(VectorMap(mapName,mapsFolder.c_str(),true));
  }
  if(maps.size()>0)
    currentMap = &maps[0];
  if(debug) printf("Done Loading Atlas.\n");
}

void VectorLocalization2D::setLocation(vector2f loc, float angle, const char* map, float locationUncertainty, float angleUncertainty)
{
  setLocation(loc, angle, locationUncertainty, angleUncertainty);

  bool found = false;
  int mapIndex=0;
  for(unsigned int i=0; !found && i<maps.size(); i++){
    if(maps[i].mapName.compare(map)==0){
      mapIndex = i;
      found = true;
    }
  }
  if(found){
    currentMap = &maps[mapIndex];
  }else{
    char buf[4096];
    snprintf(buf, 4095, "Unknown map: \"%s\"",map);
    TerminalWarning(buf);
  }  
}

void VectorLocalization2D::setLocation(vector2f loc, float angle, float locationUncertainty, float angleUncertainty)
{
  ParticleInitializer particleInitializer(loc, angle, locationUncertainty, angleUncertainty);
  PARTICLE_INITIALIZER = &particleInitializer;
  graph->transform_vertices(initializeParticle);
  PARTICLE_INITIALIZER = NULL;
}

void VectorLocalization2D::setMap(const char* map)
{
  bool found = false;
  int mapIndex=0;
  for(unsigned int i=0; !found && i<maps.size(); i++){
    if(maps[i].mapName.compare(map)==0){
      mapIndex = i;
      found = true;
    }
  }
  if(found){
    currentMap = &maps[mapIndex];
  }else{
    char buf[4096];
    snprintf(buf, 4095, "Unknown map: \"%s\"",map);
    TerminalWarning(buf);
  }
}

void VectorLocalization2D::initialize(const char* mapName, vector2f loc, float angle, float locationUncertainty, float angleUncertainty)
{
  static const bool debug = false;

  currentMap = &maps[0];
  bool found = false;
  for(unsigned int i=0; i<maps.size(); i++){
    if(maps[i].mapName.compare(mapName)==0){
      found = true;
      currentMap = &maps[i];
    }
  }
  if(!found){
    char buf[2048];
    snprintf(buf,2047,"Map %s not found in Atlas! Reverting to map %s.",mapName,maps[0].mapName.c_str());
    TerminalWarning(buf);
  }

  setLocation(loc, angle, locationUncertainty, angleUncertainty);
  computeLocation(loc, angle);

  laserEval.numCorrespondences = 0;
  laserEval.numObservedPoints = 0;
  laserEval.runTime = 0;
  laserEval.stage0Weights = 0;
  laserEval.stageRWeights = 0;

  pointCloudEval.numCorrespondences = 0;
  pointCloudEval.numObservedPoints = 0;
  pointCloudEval.runTime = 0;
  pointCloudEval.stage0Weights = 0;
  pointCloudEval.stageRWeights = 0;

  if(debug) printf("\nDone.\n");
}

void VectorLocalization2D::predict(float dx, float dy, float dtheta, const MotionModelParams &motionParams)
{
  lastDistanceMoved += vector2f(dx,dy).rotate(lastAngleTurned);
  lastAngleTurned += dtheta;

  Motion motion(dx, dy, dtheta, motionParams);
  MOTION = &motion;
  graph->transform_vertices(predictParticle);
  MOTION = NULL;
}

float VectorLocalization2D::observationWeightPointCloud(vector2f loc, float angle, const vector< vector2f >& pointCloud, const vector< vector2f >& pointNormals, const PointCloudParams & pointCloudParams) const
{
  static const bool UseAnalyticRender = true;
  static const bool debug = false;

  float logOutOfRangeProb = pointCloudParams.logOutOfRangeProb;
  float logObstacleProb = pointCloudParams.logObstacleProb;
  float logShortHitProb = pointCloudParams.logShortHitProb;
  float corelationFactor = pointCloudParams.corelationFactor;
  float minRange = pointCloudParams.minRange;
  float maxRange = pointCloudParams.maxRange;
  float sqMaxRange = sq(maxRange);
  float stdDev = pointCloudParams.stdDev;
  float curAngle;

  vector<line2f> lines;
  vector<int> lineCorrespondences;

  float a0 = angle-0.5*pointCloudParams.fieldOfView;
  float a1 = angle+0.5*pointCloudParams.fieldOfView;
  if(UseAnalyticRender){
    lineCorrespondences = currentMap->getRayToLineCorrespondences(loc, angle, a0, a1, pointCloud, minRange, maxRange, true, &lines);
  }else{
    lineCorrespondences = currentMap->getRayToLineCorrespondences(loc, angle, a0, a1, pointCloud, minRange, maxRange);
    lines = currentMap->lines;
  }

  float logTotalWeight = 0.0;

  Matrix2f rotMat1;
  rotMat1 = Rotation2Df(angle);

  int numCorrespondences = 0;
  Vector2f heading, scanPoint, lineNorm, curPoint, attraction, locE(V2COMP(loc)), rotatedNormal;

  for(int i=0; i<pointCloud.size(); i++){
    curPoint = rotMat1*Vector2f(V2COMP(pointCloud[i])) + locE;
    attraction = Vector2f(0.0,0.0);
    if(lineCorrespondences[i]>=0){
      numCorrespondences++;
      if(UseAnalyticRender){
        attraction = observationFunction(lines[lineCorrespondences[i]], curPoint);
      }else{
        attraction = observationFunction(currentMap->lines[lineCorrespondences[i]], curPoint);
      }
      //logTotalWeight += -attraction.squaredNorm()*corelationFactor/stdDev;
      /**/
      if(UseAnalyticRender){
        lineNorm = Vector2f(V2COMP(lines[lineCorrespondences[i]].Perp()));
      }else{
        lineNorm = Vector2f(V2COMP(currentMap->lines[lineCorrespondences[i]].Perp()));
      }
      rotatedNormal = rotMat1*Vector2f(V2COMP(pointNormals[i]));
      if(true || fabs(lineNorm.dot(rotatedNormal)) > pointCloudParams.minCosAngleError){
        logTotalWeight += -attraction.squaredNorm()*corelationFactor/stdDev;
      }else{
        logTotalWeight += logObstacleProb*corelationFactor;
      }
      /**/
    }else if( pointCloud[i].sqlength()<sqMaxRange){
      logTotalWeight += logObstacleProb*corelationFactor;
    }else{
      logTotalWeight += logOutOfRangeProb*corelationFactor;
    }
  }
  //return numCorrespondences;
  return exp(logTotalWeight);
}

float VectorLocalization2D::observationWeightLidar(vector2f loc, float angle, const LidarParams &lidarParams, const vector<Vector2f> &laserPoints)
{
  static const bool UseAnalyticRender = true;
  static const bool debug = false;

  float numRays = lidarParams.numRays;
  float minRange = lidarParams.minRange;
  float maxRange = lidarParams.maxRange;
  float logOutOfRangeProb = lidarParams.logOutOfRangeProb;
  float logObstacleProb = lidarParams.logObstacleProb;
  float logShortHitProb = lidarParams.logShortHitProb;
  
  vector2f laserLoc = loc + vector2f(lidarParams.laserToBaseTrans.x(),lidarParams.laserToBaseTrans.y()).rotate(angle);
  vector<line2f> lines;
  vector<int> lineCorrespondences;

  if(UseAnalyticRender){
    lineCorrespondences = currentMap->getRayToLineCorrespondences(laserLoc, angle, lidarParams.angleResolution, numRays, 0.0, maxRange, true, &lines);
  }else{
    lineCorrespondences = currentMap->getRayToLineCorrespondences(laserLoc, angle, lidarParams.angleResolution, numRays, 0.0, maxRange, false, NULL);
  }
  
  float *scanRays = lidarParams.laserScan;
  float curAngle;
  Vector2f attraction;
  float midScan = float(numRays)/2.0;
  
  Matrix2f robotAngle;
  robotAngle = Rotation2Df(angle);
  Vector2f heading, scanPoint, laserLocE(V2COMP(laserLoc));
  
  float logTotalWeight = 0.0;
  
  Vector2f scanPoint2, scanDir, lineDir;
  Matrix2f robotAngle2;
  robotAngle2 = Rotation2Df(angle+lidarParams.angleResolution);
  for(int i=0; i<numRays-1; i++){
    if(scanRays[i]>maxRange || scanRays[i]<minRange){
      logTotalWeight += logObstacleProb*lidarParams.correlationFactor;
      continue;
    }
    
    if(lineCorrespondences[i]>=0){
      scanPoint = laserLocE + robotAngle*laserPoints[i];
      scanPoint2 = laserLocE + robotAngle2*laserPoints[i+1];
      scanDir = (scanPoint2-scanPoint).normalized();
      if(UseAnalyticRender){
        lineDir = Vector2f(V2COMP(lines[lineCorrespondences[i]].Dir()));
      }else{
        lineDir = Vector2f(V2COMP(currentMap->lines[lineCorrespondences[i]].Dir()));
      }
      if(fabs(lineDir.dot(scanDir))>lidarParams.minCosAngleError){
        if(UseAnalyticRender)
          attraction = observationFunction(lines[lineCorrespondences[i]], scanPoint);
        else
          attraction = observationFunction(currentMap->lines[lineCorrespondences[i]], scanPoint);
        logTotalWeight += -min(attraction.squaredNorm()/lidarParams.lidarStdDev, -logShortHitProb)*lidarParams.correlationFactor;
      }else{
        logTotalWeight += logObstacleProb*lidarParams.correlationFactor;
      }
      
    }else if(scanRays[i]<maxRange){
      logTotalWeight += logObstacleProb*lidarParams.correlationFactor;
    }else{
      logTotalWeight += logOutOfRangeProb*lidarParams.correlationFactor;
    }
  }
  //if(debug) printf("numCorrespondences: %d\n",numCorrespondences);
  //exit(-1);
  //return numCorrespondences;
  return exp(logTotalWeight);
}

void VectorLocalization2D::updateLidar(const LidarParams &lidarParams, const MotionModelParams & motionParams)
{
  static const bool debug = false;
  static const bool usePermissibility = true;
  
  double tStart = GetTimeSec();
  
  //Transform laserpoints to robot frame
  vector< Vector2f > laserPoints(lidarParams.numRays);
  for(int i=0; i<lidarParams.numRays; i++){
    laserPoints[i] = lidarParams.laserToBaseTrans + lidarParams.laserToBaseRot*lidarParams.scanHeadings[i]*lidarParams.laserScan[i];
  }
  
  //Compute importance weights
  int N = int(numParticles);
  float totalWeight = 0.0;
  if(debug) printf("\nParticle weights:\n");
  /* FIXME -- This LIDAR update is not curently done using the graph
  for(int i=0; i<N; i++){
    Particle2D &p = particles[i];
    p.weight = observationWeightLidar(p.loc, p.angle, lidarParams, laserPoints);
    totalWeight += p.weight;
    if(debug) printf("%2d: %f\n", i, p.weight);
  }
  // Normalize weights
  for(int i=0; i<N; i++) particles[i].weight /= totalWeight;
  */

  updateTime = GetTimeSec() - tStart;
}

void VectorLocalization2D::updatePointCloud(const vector<vector2f>& pointCloud, vector<vector2f>& pointNormals, const MotionModelParams &motionParams, const PointCloudParams &pointCloudParams)
{
  static const bool debug = false;
  static const bool usePermissibility = true;

  double tStart = GetTimeSec();

  //Gather together parameters to perform update
  Update update(motionParams, pointCloud, pointNormals, pointCloudParams);
  UPDATE = &update;
  //Compute importance weights
  graph->transform_vertices(updatePointCloudParticle);
  //Compute total weight
  TOTAL_WEIGHT = graph->map_reduce_vertices<TotalWeightReducer>(TotalWeightReducer::getWeight).weight;
  //Normalize weights
  graph->transform_vertices(normalizeWeightParticle);

  updateTime = GetTimeSec() - tStart;
}

void predictParticle(graph_type::vertex_type& v)
{
  static const bool debug = false;

  if(debug)
  {
    printf("predict before: %7.2f,%7.2f %6.2f\u00b0 ", v.data().loc.x, v.data().loc.y,
        DEG(v.data().angle));
  }

  //Motion model
  vector2f delta(MOTION->dx, MOTION->dy);
  float d_trans = delta.length();

  //Predict rotation
  float dtheta = MOTION->dtheta;
  dtheta += randn(MOTION->motionParams->Alpha1*dtheta, 0.0f) + 
            randn(MOTION->motionParams->Alpha2*d_trans, 0.0f);
  v.data().angle = angle_mod(v.data().angle + dtheta);

  //Predict translation
  if(d_trans>FLT_MIN){
    delta = delta.rotate(v.data().angle);
    delta = delta.norm(d_trans + randn(MOTION->motionParams->Alpha3*d_trans, 0.0f));
  }

  if(debug)
    printf("delta: %f,%f %f\u00b0 ", V2COMP(delta), DEG(dtheta));

  v.data().loc.x += delta.x;
  v.data().loc.y += delta.y;

  if(debug)
    printf("after: %7.2f,%7.2f %6.2f\u00b0\n", v.data().loc.x, v.data().loc.y, DEG(v.data().angle));
}

inline Vector2f VectorLocalization2D::attractorFunction(line2f l, Vector2f p, float attractorRange, float margin) const
{
  static const bool debug = false;
  
  Vector2f attraction(0.0,0.0), dir(V2COMP(l.Dir())), p0(V2COMP(l.P0())), p1(V2COMP(l.P1()));
  
  float location = (p-p0).dot(dir);
  
  if(location<-margin || location>l.Length()+margin){
    return attraction;
  }
  
  attraction = (p-p0) - dir*location ;
  
  float d = attraction.norm();
  /*
  if(d>0.5*attractorRange){
    float d2 = max(0.0f, attractorRange - d);
    attraction *= 2.0f*d2/attractorRange;
  }
  */
  if(d>attractorRange){
    attraction = Vector2f::Zero();
  }
  return attraction;
}

inline Vector2f VectorLocalization2D::observationFunction(line2f l, Vector2f p) const
{
  static const bool debug = false;
  Vector2f attraction(0.0,0.0), dir(V2COMP(l.Dir())), p0(V2COMP(l.P0())), p1(V2COMP(l.P1()));
  float location = (p-p0).dot(dir);
  attraction = (p-p0) - dir*location ;
  return attraction;
}

void VectorLocalization2D::getPointCloudGradient(int particleIdx, vector2f loc, float angle, vector2f& locGrad, float& angleGrad, const std::vector< vector2f >& pointCloud, const std::vector< vector2f >& pointNormals, float& logWeight, const VectorLocalization2D::PointCloudParams& pointCloudParams, const vector<int> & lineCorrespondences, const vector<line2f> &lines) const
{
  static const bool UseAnalyticRender = false;
  static const bool debug = false;
  static const bool EnableProfiling = false;

  FunctionTimer *ft;
  if(EnableProfiling)
    ft = new FunctionTimer(__FUNCTION__);
  
  float numTotalPoints = pointCloud.size();
  
  if(EnableProfiling) ft->Lap(__LINE__);
  //vector<int> lineCorrespondences;
  //vector<line2f> lines;
  
  /*
  float minRange = pointCloudParams.minRange;
  float maxRange = pointCloudParams.maxRange;
  float a0 = angle-0.5*pointCloudParams.fieldOfView;
  float a1 = angle+0.5*pointCloudParams.fieldOfView;
  if(UseAnalyticRender){
    lineCorrespondences = currentMap->getRayToLineCorrespondences(loc, angle, a0, a1, pointCloud, minRange, maxRange, true, &lines);
  }else{
    lineCorrespondences = currentMap->getRayToLineCorrespondences(loc, angle, a0, a1, pointCloud, minRange, maxRange);
    lines = currentMap->lines;
  }
  */
  if(EnableProfiling) ft->Lap(__LINE__);
  
  float numPoints = 0;
  
  float cosAngle;
  int noCorrespondences = 0;
  logWeight = 0.0;
  
  int numObservedPoints = int(pointCloud.size());
  
  //Construct gradients per point in point cloud

  std::vector<Vector2f> gradients2(numObservedPoints);
  std::vector<Vector2f> points2(numObservedPoints);
  int numPointsInt = 0;
  
  Vector2f curPoint, locE(V2COMP(loc)), attraction, lineNorm, rotatedNormal;
  Matrix2f rotMat1;
  rotMat1 = Rotation2Df(angle);
  
  for(int i=0; i<numObservedPoints; i++){
    curPoint = rotMat1*Vector2f(V2COMP(pointCloud[i])) + locE;
    rotatedNormal = rotMat1*Vector2f(V2COMP(pointNormals[i]));
    attraction = Vector2f(0,0);

    if(lineCorrespondences[i]>=0 && lineCorrespondences[i]<lines.size()){
      lineNorm = Vector2f(V2COMP(lines[lineCorrespondences[i]].Perp()));
      cosAngle = fabs(lineNorm.dot(rotatedNormal));

      //Attraction only for valid correspondences
      if(cosAngle > pointCloudParams.minCosAngleError){
        attraction = attractorFunction(lines[lineCorrespondences[i]], curPoint,pointCloudParams.attractorRange, pointCloudParams.correspondenceMargin);
        //Add the point and attraction (only if non-zero)
        //logWeight += -min(attraction.squaredNorm()/pointCloudParams.stdDev, -pointCloudParams.logShortHitProb)*pointCloudParams.corelationFactor;
        gradients2[numPointsInt] = attraction;
        points2[numPointsInt] = curPoint;
        numPointsInt++;
      }
    }else{
      noCorrespondences++;
    }
  }
  gradients2.resize(numPointsInt);
  points2.resize(numPointsInt);
  numPoints = float(numPointsInt);

  if(debug) printf("No correspondences: %d/%d \n",noCorrespondences,int(pointCloud.size()));

  if(numPoints<pointCloudParams.minPoints){
    locGrad.zero();
    angleGrad = 0.0;
    return;
  }

  //Estimate translation and rotation
  Vector2f locGradE(0,0);
  Vector2f heading(0,0), curHeading, r;
  float headingAngle;
  for(int i = 0; i<numPointsInt; i++){
    r = points2[i] - locE;
    if(r.squaredNorm()<sq(0.001))
      continue;

    locGradE += gradients2[i];
    headingAngle = eigenCross(r, gradients2[i]);
    curHeading = Vector2f(cos(headingAngle),sin(headingAngle));
    heading += curHeading;
  }
  locGradE = locGradE/numPoints;
  locGrad.set(locGradE.x(),locGradE.y());
  heading = heading/numPoints;

  angleGrad = bound(atan2(heading.y(),heading.x()),-pointCloudParams.maxAngleGradient,pointCloudParams.maxAngleGradient);

  locGrad = locGrad.bound(pointCloudParams.maxLocGradient);
  if(debug) printf("LocGrad: %6.2f %6.2f AngleGrad:%6.1f\u00b0\n",V2COMP(locGrad), DEG(angleGrad));

  if(EnableProfiling) delete ft;
}

void VectorLocalization2D::getLidarGradient(vector2f loc, float angle, vector2f& locGrad, float& angleGrad, float& logWeight, VectorLocalization2D::LidarParams lidarParams, const vector<Vector2f>& laserPoints, const vector<int> & lineCorrespondences, const vector<line2f> &lines)
{
  static const bool debug = false;
  static const bool EnableProfiling = false;
  static const bool UseAnalyticRender = true;
  
  FunctionTimer *ft;
  if(EnableProfiling)
    ft = new FunctionTimer(__PRETTY_FUNCTION__);
  
  float logShortHitProb = lidarParams.logShortHitProb;
  Matrix2f robotAngle;
  robotAngle = Rotation2Df(angle);
  Vector2f laserLocE = Vector2f(V2COMP(loc)) + robotAngle*(lidarParams.laserToBaseTrans);
  vector2f laserLoc(laserLocE.x(), laserLocE.y());
  
  if(EnableProfiling) ft->Lap(__LINE__);
  
  /*
  if(UseAnalyticRender){
    lineCorrespondences = currentMap->getRayToLineCorrespondences(laserLoc, angle, lidarParams.angleResolution, lidarParams.numRays, lidarParams.minRange, lidarParams.maxRange, true, &lines);
  }else{
    lineCorrespondences = currentMap->getRayToLineCorrespondences(laserLoc, angle, lidarParams.angleResolution, lidarParams.numRays, lidarParams.minRange, lidarParams.maxRange);
  }
  */
  if(EnableProfiling) ft->Lap(__LINE__);
  
  points.clear();
  gradients.clear();

  float *scanRays = lidarParams.laserScan;
  float minRange = lidarParams.minRange;
  float maxRange = lidarParams.maxRange;
  float curAngle;
  Vector2f heading, heading2;
  float numPoints = 0;
  float midScan = float(lidarParams.numRays)/2.0;
  float cosAngle;
  int noCorrespondences = 0;
  logWeight = 0.0;
  curAngle = angle - midScan*lidarParams.angleResolution;
  //Construct gradients per point in point cloud
  if(EnableProfiling) ft->Lap(__LINE__);
  
  Matrix2f robotAngle2;
  robotAngle2 = Rotation2Df(angle+lidarParams.angleResolution);
  Vector2f scanPoint, scanPoint2, scanDir, lineDir, attraction;
  
  //TODO: Speed up this section!!
  //===========================================================================
  laserEval.numObservedPoints = 0;
  for(int i=0; i<lidarParams.numRays-1; i++, curAngle+=lidarParams.angleResolution){
    if(scanRays[i]<minRange || scanRays[i]>maxRange)
      continue;
    laserEval.numObservedPoints++;  
    scanPoint = laserLocE + robotAngle*laserPoints[i];
    
    if(lineCorrespondences[i]>=0){
      scanPoint2 = laserLocE + robotAngle2*laserPoints[i+1];
      scanDir = (scanPoint2-scanPoint).normalized();
      if(UseAnalyticRender){
        lineDir = Vector2f(V2COMP(lines[lineCorrespondences[i]].Dir()));
        attraction = attractorFunction(lines[lineCorrespondences[i]], scanPoint, lidarParams.attractorRange, lidarParams.correspondenceMargin);
      }else{
        lineDir = Vector2f(V2COMP(currentMap->lines[lineCorrespondences[i]].Dir()));
        attraction = attractorFunction(currentMap->lines[lineCorrespondences[i]], scanPoint, lidarParams.attractorRange, lidarParams.correspondenceMargin);
      }
      logWeight += -min(attraction.squaredNorm()/lidarParams.lidarStdDev, -logShortHitProb)*lidarParams.correlationFactor;
      cosAngle = fabs(lineDir.dot(scanDir));
      if(cosAngle>lidarParams.minCosAngleError){
        gradients.push_back(attraction);
        points.push_back(scanPoint);
      }
    }else{
      noCorrespondences++;
    }
    
  }  
  numPoints = float(points.size());
  laserEval.numCorrespondences = int(points.size());
  
  //if(debug) printf("No correspondences: %d/%d numPoints:%d\n",noCorrespondences,lidarParams.numRays,int(numPoints));
  //===========================================================================
  
  if(EnableProfiling) ft->Lap(__LINE__);
  
  if(numPoints<lidarParams.minPoints){
    locGrad.zero();
    angleGrad = 0.0;
    if(EnableProfiling) delete ft;
    return;
  }
  
  //Estimate translation and rotation
  locGrad.zero();
  heading = Vector2f(0,0);
  float headingAngle;
  laserEval.meanSqError = 0.0;
  Vector2f locE(V2COMP(loc)), r, locGradE(0.0,0.0);
  for(int i = 0; i<gradients.size(); i++){
    r = points[i]-locE;
    laserEval.meanSqError += gradients[i].squaredNorm();
    if(r.squaredNorm()<sq(0.03))
      continue;
    locGradE += gradients[i];
    headingAngle = eigenCross(r, gradients[i]);
    heading += Vector2f(cos(headingAngle),sin(headingAngle));
  }
  locGradE /= numPoints;
  locGrad.set(locGradE.x(),locGradE.y());
  heading /= numPoints;
  laserEval.meanSqError /= numPoints;
  
  if(EnableProfiling) ft->Lap(__LINE__);
    
  locGrad = locGrad.bound(lidarParams.maxLocGradient);
  angleGrad = bound(atan2(heading.y(),heading.x()),-lidarParams.maxAngleGradient,lidarParams.maxAngleGradient);
  if(EnableProfiling) delete ft;
  if(debug) printf("Gradient: %.4f,%.4f %.2f\u00b0\n",V2COMP(locGrad),DEG(angleGrad));
}

void VectorLocalization2D::refineLocationLidar(vector2f& loc, float& angle, float& initialWeight, float& finalWeight, const LidarParams &lidarParams, const vector<Vector2f> &laserPoints)
{
  static const bool debug = false;
  static const bool UseAnalyticRender = true;
  
  //Do gradient descent for this particle
  vector2f locGrad(0.0,0.0);
  float angleGrad = 0.0;
  bool beingRefined = true;
  float weight;
  
  if(debug) printf("before: %.4f,%.4f %.2f\u00b0\n",V2COMP(loc),DEG(angle));
  
  Matrix2f robotAngle;
  robotAngle = Rotation2Df(angle);  
  Vector2f laserLocE = Vector2f(V2COMP(loc)) + robotAngle*(lidarParams.laserToBaseTrans);
  vector2f laserLoc(laserLocE.x(), laserLocE.y());
  vector<line2f> lines;
  vector<int> lineCorrespondences;

  if(UseAnalyticRender){
    lineCorrespondences = currentMap->getRayToLineCorrespondences(laserLoc, angle, lidarParams.angleResolution, lidarParams.numRays, lidarParams.minRange, lidarParams.maxRange, true, &lines);
  }else{
    lineCorrespondences = currentMap->getRayToLineCorrespondences(laserLoc, angle, lidarParams.angleResolution, lidarParams.numRays, lidarParams.minRange, lidarParams.maxRange);
  }
  
  for(int i=0; beingRefined && i<lidarParams.numSteps; i++){
    getLidarGradient(loc,angle,locGrad,angleGrad,weight,lidarParams, laserPoints, lineCorrespondences, lines);
    if(i==0) initialWeight = exp(weight);
    loc -= lidarParams.etaLoc*locGrad;
    angle -= lidarParams.etaAngle*angleGrad;
    beingRefined = fabs(angleGrad)>lidarParams.minRefineFraction*lidarParams.maxAngleGradient && locGrad.sqlength()>sq(lidarParams.minRefineFraction*lidarParams.maxLocGradient);
  }
  
  if(debug) printf("after: %.4f,%.4f %.2f\u00b0\n",V2COMP(loc),DEG(angle));
  finalWeight = exp(weight);
}

void VectorLocalization2D::refineLocationPointCloud(int particleIdx, vector2f& loc, float& angle, const std::vector< vector2f >& pointCloud, const std::vector< vector2f >& pointNormals, const VectorLocalization2D::PointCloudParams& pointCloudParams) const
{
  static const bool debug = false;
  static const bool UseAnalyticRender = true;

  //Do gradient descent for this particle
  vector2f locGrad(0.0,0.0);
  float angleGrad = 0.0;
  bool beingRefined = true;

  float a0 = angle-0.5*pointCloudParams.fieldOfView;
  float a1 = angle+0.5*pointCloudParams.fieldOfView;
  vector<line2f> lines;
  vector<int> lineCorrespondences;

  if(UseAnalyticRender){
    lineCorrespondences = currentMap->getRayToLineCorrespondences(loc, angle, a0, a1, pointCloud,
        pointCloudParams.minRange, pointCloudParams.maxRange, true, &lines);
  }else{
    lineCorrespondences = currentMap->getRayToLineCorrespondences(loc, angle, a0, a1, pointCloud,
        pointCloudParams.minRange, pointCloudParams.maxRange);
    lines = currentMap->lines;
  }

  vector2f locPrev = loc;
  float anglePrev = angle;
  float weight;
  for(int i=0; beingRefined && i<pointCloudParams.numSteps; i++){
    getPointCloudGradient(particleIdx, loc, angle, locGrad, angleGrad, pointCloud, pointNormals, weight, pointCloudParams, lineCorrespondences, lines);
    loc -= pointCloudParams.etaLoc*locGrad;
    angle -= pointCloudParams.etaAngle*angleGrad;
    beingRefined = fabs(angleGrad)>pointCloudParams.minRefineFraction*pointCloudParams.maxAngleGradient && locGrad.sqlength()>sq(pointCloudParams.minRefineFraction*pointCloudParams.maxLocGradient);
  }
  if(debug) printf("gradient: %7.3f,%7.3f %6.1f\u00b0\n",V2COMP(loc-locPrev),DEG(angle-anglePrev));
}

void VectorLocalization2D::refineLidar(const LidarParams &lidarParams)
{
  double tStart = GetTimeSec();
  laserEval.stage0Weights = 0.0;
  laserEval.stageRWeights = 0.0;
  laserEval.lastRunTime = GetTimeSec();
  
  // Transform laserpoints to robot frame
  vector< Vector2f > laserPoints(lidarParams.numRays);
  for(int i=0; i<lidarParams.numRays; i++){
    laserPoints[i] = lidarParams.laserToBaseTrans + lidarParams.laserToBaseRot*lidarParams.scanHeadings[i]*lidarParams.laserScan[i];
  }
  
  /* FIXME -- The LIDAR update is not currently done using the graph
  if(lidarParams.numSteps>0){
    for(int i=0; i<numParticles; i++){
      float initialWeight, finalWeight;
      refineLocationLidar(particles[i].loc, particles[i].angle, initialWeight, finalWeight, lidarParams, laserPoints);
      laserEval.stage0Weights += initialWeight;
      laserEval.stageRWeights += finalWeight;
    }
  }
  */
  refineTime = GetTimeSec() - tStart;
  laserEval.runTime = refineTime;
}


void VectorLocalization2D::reducePointCloud(const std::vector< vector2f >& pointCloud, const std::vector< vector2f >& pointNormals, vector< vector2f >& reducedPointCloud, vector< vector2f >& reducedPointNormals)
{
  static const bool debug = false;
  static const float eps = -RAD(0.001);

  vector<pair<float, int> > angles;
  size_t N = pointCloud.size();
  pair<float,int> valuePair;

  for(unsigned int i=0; i<N; i++){
    valuePair.first = pointCloud[i].angle();
    valuePair.second = i;
    angles.push_back(valuePair);
  }
  //N = min(angles.size(),N);
  sort<vector<pair<float, int> >::iterator >(angles.begin(), angles.end());

  reducedPointCloud.resize(N);
  reducedPointNormals.resize(N);
  int j=0;
  float lastAngle;
  for(unsigned int i=0; i<N; ){
    reducedPointCloud[j] = pointCloud[angles[i].second];
    reducedPointNormals[j] = pointNormals[angles[i].second];
    j++;
    lastAngle = angles[i].first;

    for(i++;i<N && angles[i].first<=lastAngle+eps;){
      i++;
    }
  }
  N = j;
  reducedPointCloud.resize(N);
  reducedPointNormals.resize(N);

  if(debug){
    printf("\n");
    for(unsigned int i=0; i<N; i++){
      //printf("%4d: %6.2f\u00b0 %4d\n",int(i),DEG((*it).first),(*it).second);
      printf("%4d: %6.2f\u00b0\n",int(i),DEG(angles[i].first));
    }
  }
}

void VectorLocalization2D::refinePointCloud(const vector<vector2f> &pointCloud, const vector<vector2f> &pointNormals, const PointCloudParams &pointCloudParams)
{
  //FunctionTimer ft(__PRETTY_FUNCTION__);
  double tStart = GetTimeSec();

  if(pointCloudParams.numSteps>0){
    Refinement refinement(pointCloud, pointNormals, pointCloudParams);
    REFINEMENT = &refinement;
    graph->transform_vertices(refinePointCloudParticle);
    REFINEMENT = NULL;
  }

  refineTime = GetTimeSec() - tStart;
  pointCloudEval.runTime = refineTime;
}

void refinePointCloudParticle(graph_type::vertex_type& v)
{
  localization->refineLocationPointCloud(v.id(), v.data().loc, v.data().angle, *REFINEMENT->pointCloud,
      *REFINEMENT->pointNormals, *REFINEMENT->pointCloudParams);
}

void updatePointCloudParticle(graph_type::vertex_type& v)
{
  static const bool debug = false;

  if(debug) printf("\nParticle weights:\n");
  //Compute importance weights
  v.data().weight = localization->observationWeightPointCloud(v.data().loc, v.data().angle, *UPDATE->pointCloud, *UPDATE->pointNormals, *UPDATE->pointCloudParams);
  if(debug) printf("%2lu: %f\n", v.id(), v.data().weight);
}

void normalizeWeightParticle(graph_type::vertex_type& v) 
{
  v.data().weight /= TOTAL_WEIGHT;
}

void VectorLocalization2D::computeLocation(vector2f& loc, float& angle)
{
  //sum all the poses (taking care od adding up headings appropriately)
  PoseReducer r = graph->map_reduce_vertices<PoseReducer>(PoseReducer::getPose);
  //maximum likelihood estimation of the pose
  loc.x = r.x/numParticles;
  loc.y = r.y/numParticles;
  angle = std::atan2(r.heading_y, r.heading_x);
}

void VectorLocalization2D::resample(Resample type)
{
  switch(type){
  case MultinomialResampling:
    multinomialResample();
    break;
  }
}

void VectorLocalization2D::multinomialResample()
{
  engine->signal_all();
  engine->start();
}

void VectorLocalization2D::saveProfilingStats(FILE* f)
{
  fprintf(f, "%f, %f, ",refineTime, updateTime);
}

void VectorLocalization2D::getEvalValues(VectorLocalization2D::EvalValues& _laserEval, VectorLocalization2D::EvalValues& _pointCloudEval)
{
  laserEval.stage0Weights /= float(numParticles);
  laserEval.stageRWeights /= float(numParticles);
  pointCloudEval.stage0Weights /= float(numParticles);
  pointCloudEval.stageRWeights /= float(numParticles);
  _laserEval = laserEval;
  _pointCloudEval = pointCloudEval;
}

void VectorLocalization2D::getUncertainty(float& _angleUnc, float& _locUnc)
{
  _angleUnc = currentAngleStdDev;
  _locUnc = currentLocStdDev.length();
}

VectorLocalization2D::Motion::Motion()
{
  dx = 0.0;
  dy = 0.0;
  dtheta = 0.0;
  motionParams = NULL;
}

VectorLocalization2D::Motion::Motion(float _dx, float _dy, float _dtheta, const MotionModelParams &_motionParams)
{
  dx = _dx;
  dy = _dy;
  dtheta = _dtheta;
  motionParams = &_motionParams;
}

VectorLocalization2D::Refinement::Refinement()
{
  pointCloud = NULL;
  pointNormals = NULL;
  pointCloudParams = NULL;
}

VectorLocalization2D::Refinement::Refinement(const std::vector<vector2f>& _pointCloud,
    const std::vector<vector2f>& _pointNormals, const PointCloudParams& _pointCloudParams)
{
  pointCloud = &_pointCloud;
  pointNormals = &_pointNormals;
  pointCloudParams = &_pointCloudParams;
}

VectorLocalization2D::Update::Update()
{
  motionParams = NULL;
  pointCloud = NULL;
  pointNormals = NULL;
  pointCloudParams = NULL;
}

VectorLocalization2D::Update::Update(const MotionModelParams& _motionParams, const std::vector<vector2f>& _pointCloud,
    const std::vector<vector2f>& _pointNormals, const PointCloudParams& _pointCloudParams)
{
  motionParams = &_motionParams;
  pointCloud = &_pointCloud;
  pointNormals = &_pointNormals;
  pointCloudParams = &_pointCloudParams;
}

PoseReducer PoseReducer::getPose(const graph_type::vertex_type& v) {
  PoseReducer r;

  r.x = v.data().loc.x;
  r.y = v.data().loc.y;

  r.heading_x = std::cos(v.data().angle);
  r.heading_y = std::sin(v.data().angle);

  return r;
}

PoseReducer& PoseReducer::operator+=(const PoseReducer& other) {
  x += other.x;
  y += other.y;

  heading_x += other.heading_x;
  heading_y += other.heading_y;

  return *this;
}

TotalWeightReducer TotalWeightReducer::getWeight(const graph_type::vertex_type& v) {
  return TotalWeightReducer(v.data().weight);
}

TotalWeightReducer& TotalWeightReducer::operator+=(const TotalWeightReducer& other) {
  weight += other.weight;
  return *this;
}

void initializeParticle(graph_type::vertex_type& v) {
  v.data() = Particle2D(
      randn(PARTICLE_INITIALIZER->locationUncertainty, PARTICLE_INITIALIZER->loc.x),
      randn(PARTICLE_INITIALIZER->locationUncertainty, PARTICLE_INITIALIZER->loc.y),
      randn(PARTICLE_INITIALIZER->angleUncertainty, PARTICLE_INITIALIZER->angle),
      1.0);
}

VectorLocalization2D::ParticleInitializer::ParticleInitializer() {
  loc.zero();
  angle = 0.0;
  locationUncertainty = 0.0;
  angleUncertainty = 0.0;
}

VectorLocalization2D::ParticleInitializer::ParticleInitializer(vector2f _loc, float _angle, float _locationUncertainty, float _angleUncertainty) {
  loc.x = _loc.x;
  loc.y = _loc.y;
  angle = _angle;
  locationUncertainty = _locationUncertainty;
  angleUncertainty = _angleUncertainty;
}

Resampler::Resampler() { }

Resampler::Resampler(const Particle2D& particle) {
  in_particles.push_back(particle);
}

Resampler& Resampler::operator+=(const Resampler& other) {
  for (unsigned int i = 0; i < other.in_particles.size(); i++)
    in_particles.push_back(other.in_particles[i]);

  return *this;
}

void Resampler::save(graphlab::oarchive& oarc) const {
  oarc << in_particles;
}

void Resampler::load(graphlab::iarchive& iarc) {
  iarc >> in_particles;
}

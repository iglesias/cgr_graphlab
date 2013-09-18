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

static const bool UseAnalyticRender = true;

std::vector<line2f> debugLines;

//Pointer to the parameters required to initialize particles
const VectorLocalization2D::ParticleInitializer* PARTICLE_INITIALIZER = NULL;
//Pointer to the motion currently processed
const VectorLocalization2D::Motion* MOTION = NULL;
//Pointer to the current refinement parameters
const VectorLocalization2D::Refinement* REFINEMENT = NULL;
//Pointer to the curent update parameters
const VectorLocalization2D::Update* UPDATE = NULL;

//Set of particles set during graph transformations
std::vector<Particle2D> particles;
//Densities computed for each particle to perform update
std::vector<float> SAMPLING_DENSITY;
//Total sampling density used for normalization
float TOTAL_DENSITY;

//Per-particle statistics of performance
std::vector<VectorLocalization2D::EvalValues> PARTICLE_POINT_CLOUD_EVAL;

//Per-particle debugging variables for the GUI
std::vector< std::vector<Vector2f> > GRADIENTS2;
std::vector< std::vector<Vector2f> > POINTS2;
std::vector<vector2f> locCorrectionP0;
std::vector<vector2f> locCorrectionP1;

//This guy is declared and created in localization_main
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

VectorLocalization2D::VectorLocalization2D(int _numParticles, graph_type& _graph, const char* _mapsFolder)
{
  assert(_numParticles > 0);

  mapsFolder = string(_mapsFolder);
  loadAtlas();
  numParticles = _numParticles;
  particles.resize(_numParticles);
  SAMPLING_DENSITY.resize(_numParticles);

  //create distributed graph
  for (int id = 0; id < numParticles; ++id)
    _graph.add_vertex(id, Particle2D());

  //commit the distributed graph, denoting that it is no longer to be modified
  _graph.finalize();

  graph = &_graph;
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
  for(unsigned int i=0; i<particles.size(); i++){
    particles[i].loc = vector2f(randn(locationUncertainty, loc.x), randn(locationUncertainty, loc.y));
    particles[i].angle = randn(angleUncertainty, angle);
  }

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
  for(unsigned int i=0; i<particles.size(); i++){
    particles[i].loc = vector2f(randn(locationUncertainty, loc.x), randn(locationUncertainty, loc.y));
    particles[i].angle = randn(angleUncertainty, angle);
  }
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
  
  ParticleInitializer particleInitializer(loc, angle, locationUncertainty, angleUncertainty);
  PARTICLE_INITIALIZER = &particleInitializer;
  graph->transform_vertices(initializeParticle);
  PARTICLE_INITIALIZER = NULL;
  
  computeLocation(loc, angle);
  
  laserEval.numCorrespondences = 0;
  laserEval.numObservedPoints = 0;
  laserEval.runTime = 0;
  laserEval.stage0Weights = 0;
  laserEval.stageRWeights = 0;
  
  PARTICLE_POINT_CLOUD_EVAL.resize(numParticles);
  pointCloudEval.numCorrespondences = 0;
  pointCloudEval.numObservedPoints = 0;
  pointCloudEval.runTime = 0;
  pointCloudEval.stage0Weights = 0;
  pointCloudEval.stageRWeights = 0;

  GRADIENTS2.resize(numParticles);
  POINTS2.resize(numParticles);
  locCorrectionP0.resize(numParticles);
  locCorrectionP1.resize(numParticles);
  
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

float VectorLocalization2D::motionModelWeight(vector2f loc, float angle, const MotionModelParams &motionParams) const
{
  float sqDensityKernelSize = sq(motionParams.kernelSize);
  float w = 0.0;
  float incr = 1.0/float(numParticles);
  for(int i=0; i<numParticles; i++){
    if( (loc-particles[i].loc).sqlength() < sqDensityKernelSize )
      w += incr;
  }
  return w;
}

float VectorLocalization2D::observationWeightPointCloud(vector2f loc, float angle, const vector< vector2f >& pointCloud, const vector< vector2f >& pointNormals, const PointCloudParams & pointCloudParams) const
{
  //static const bool UseAnalyticRender = false;
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
  
  // Transform laserpoints to robot frame
  vector< Vector2f > laserPoints(lidarParams.numRays);
  for(int i=0; i<lidarParams.numRays; i++){
    laserPoints[i] = lidarParams.laserToBaseTrans + lidarParams.laserToBaseRot*lidarParams.scanHeadings[i]*lidarParams.laserScan[i];
  }
  
  //Compute the sampling density
  float sqDensityKernelSize = sq(lidarParams.kernelSize);
  float totalDensity = 0.0;
  int N = int(particles.size());
  static vector<float> samplingDensity;
  if(samplingDensity.size()!=N)
    samplingDensity.resize(N);
  if(debug) printf("\nParticle samplingDensity:\n");
  for(int i=0; i<N; i++){
    float w = 0.99;
    for(int j=0; j<N; j++){
      if(i==j)
        continue;
      if( (particles[j].loc - particles[i].loc).sqlength() < sqDensityKernelSize && fabs(angle_diff(particles[j].angle, particles[i].angle))<RAD(20.0))
        w++;
    }
    samplingDensity[i] = w;
    totalDensity += w;
    if(debug) printf("%2d:%f\n",i,w);
  }
  // Normalize densities, not really necessary since resampling does not require normalized weights
  for(int i=0; i<N; i++){
    samplingDensity[i] /= totalDensity;
  }
  
  //Compute importance weights = observation x motion / samplingDensity
  if(debug) printf("\nParticle weights:\n");
  for(int i=0; i<N; i++){
    Particle2D &p = particles[i];
    float w1 = observationWeightLidar(p.loc, p.angle, lidarParams, laserPoints);
    float w2 = motionModelWeight(p.loc, p.angle, motionParams);
    p.weight = w1*w2/samplingDensity[i];
    if(debug) printf("%2d: %f , %f , %f\n",i,w1,w2,p.weight);
  }
  
  updateTime = GetTimeSec() - tStart;
}

void VectorLocalization2D::updatePointCloud(const vector<vector2f>& pointCloud, vector<vector2f>& pointNormals, const MotionModelParams &motionParams, const PointCloudParams &pointCloudParams)
{
  static const bool debug = false;
  static const bool usePermissibility = true;

  double tStart = GetTimeSec();

  Update update(motionParams, pointCloud, pointNormals, pointCloudParams);
  UPDATE = &update;

  graph->transform_vertices(samplingDensityPointCloudParticle);

  SamplingDensityReducer r = graph->map_reduce_vertices<SamplingDensityReducer>(SamplingDensityReducer::computeTotalDensity);
  // Compute total density to later normalize densities, not really necessary though since resampling
  // does not require normalized weights
  TOTAL_DENSITY = r.samplingDensity;

  graph->transform_vertices(updatePointCloudParticle);

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

  particles[v.id()] = v.data();

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
  if(debug){
    debugLines.push_back( line2f( vector2f(p.x(),p.y()) ,vector2f((p-attraction).x(),(p-attraction).y()) ) );
  }
  return attraction;
}

inline Vector2f VectorLocalization2D::observationFunction(line2f l, Vector2f p) const
{
  static const bool debug = false;
  Vector2f attraction(0.0,0.0), dir(V2COMP(l.Dir())), p0(V2COMP(l.P0())), p1(V2COMP(l.P1()));
  float location = (p-p0).dot(dir);
  attraction = (p-p0) - dir*location ;
  if(debug){
    debugLines.push_back( line2f( vector2f(p.x(),p.y()) ,vector2f((p-attraction).x(),(p-attraction).y()) ) );
    const float crossSize = 0.002;
    debugLines.push_back( line2f( vector2f(p.x()+crossSize,p.y()) , vector2f(p.x()-crossSize,p.y())) );
    debugLines.push_back( line2f( vector2f(p.x(),p.y()+crossSize) , vector2f(p.x(),p.y()-crossSize)) );
  }
  return attraction;
}

void VectorLocalization2D::getPointCloudGradient(int particleIdx, vector2f loc, float angle, vector2f& locGrad, float& angleGrad, const std::vector< vector2f >& pointCloud, const std::vector< vector2f >& pointNormals, float& logWeight, const VectorLocalization2D::PointCloudParams& pointCloudParams, const vector<int> & lineCorrespondences, const vector<line2f> &lines) const
{
  //static const bool UseAnalyticRender = false;
  static const bool debug = false;
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
  assert(PARTICLE_POINT_CLOUD_EVAL.size() > particleIdx);
  PARTICLE_POINT_CLOUD_EVAL[particleIdx].numObservedPoints = numObservedPoints;
  
  //Construct gradients per point in point cloud
  assert(GRADIENTS2.size() > particleIdx);
  assert(POINTS2.size() > particleIdx);

  GRADIENTS2[particleIdx].resize(numObservedPoints);
  POINTS2[particleIdx].resize(numObservedPoints);
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
        
        GRADIENTS2[particleIdx][numPointsInt] = attraction;
        POINTS2[particleIdx][numPointsInt] = curPoint;
        numPointsInt++;
      }
    }else{
      noCorrespondences++;
    }
  }
  GRADIENTS2[particleIdx].resize(numPointsInt);
  POINTS2[particleIdx].resize(numPointsInt);
  numPoints = float(numPointsInt);
  PARTICLE_POINT_CLOUD_EVAL[particleIdx].numCorrespondences = int(POINTS2[particleIdx].size());
  
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
  PARTICLE_POINT_CLOUD_EVAL[particleIdx].meanSqError = 0.0;
  for(int i = 0; i<numPointsInt; i++){
    r = POINTS2[particleIdx][i] - locE;
    PARTICLE_POINT_CLOUD_EVAL[particleIdx].meanSqError += GRADIENTS2[particleIdx][i].squaredNorm();
    if(r.squaredNorm()<sq(0.001))
      continue;
    
    locGradE += GRADIENTS2[particleIdx][i];
    headingAngle = eigenCross(r, GRADIENTS2[particleIdx][i]);
    curHeading = Vector2f(cos(headingAngle),sin(headingAngle));
    heading += curHeading;
  }
  locGradE = locGradE/numPoints;
  locGrad.set(locGradE.x(),locGradE.y());
  heading = heading/numPoints;
  PARTICLE_POINT_CLOUD_EVAL[particleIdx].meanSqError = PARTICLE_POINT_CLOUD_EVAL[particleIdx].meanSqError/numPoints;
  
  angleGrad = bound(atan2(heading.y(),heading.x()),-pointCloudParams.maxAngleGradient,pointCloudParams.maxAngleGradient);
  
  locGrad = locGrad.bound(pointCloudParams.maxLocGradient);
  if(debug) printf("LocGrad: %6.2f %6.2f AngleGrad:%6.1f\u00b0\n",V2COMP(locGrad), DEG(angleGrad));
    
  if(EnableProfiling) delete ft;
}

void VectorLocalization2D::getLidarGradient(vector2f loc, float angle, vector2f& locGrad, float& angleGrad, float& logWeight, VectorLocalization2D::LidarParams lidarParams, const vector<Vector2f>& laserPoints, const vector<int> & lineCorrespondences, const vector<line2f> &lines)
{
  static const bool EnableProfiling = false;
  
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
  static const bool debug = false;
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

void VectorLocalization2D::refineLocationPointCloud(int particleIdx, vector2f& loc, float& angle, float& initialWeight, float& finalWeight, const std::vector< vector2f >& pointCloud, const std::vector< vector2f >& pointNormals, const VectorLocalization2D::PointCloudParams& pointCloudParams) const
{
  static const bool debug = false;
  //FunctionTimer ft(__PRETTY_FUNCTION__);
  
  //Do gradient descent for this particle
  vector2f locGrad(0.0,0.0);
  float angleGrad = 0.0;
  bool beingRefined = true;
  
  float a0 = angle-0.5*pointCloudParams.fieldOfView;
  float a1 = angle+0.5*pointCloudParams.fieldOfView;
  vector<line2f> lines;
  vector<int> lineCorrespondences;

  if(UseAnalyticRender){
    lineCorrespondences = currentMap->getRayToLineCorrespondences(loc, angle, a0, a1, pointCloud, pointCloudParams.minRange, pointCloudParams.maxRange, true, &lines);
  }else{
    lineCorrespondences = currentMap->getRayToLineCorrespondences(loc, angle, a0, a1, pointCloud, pointCloudParams.minRange, pointCloudParams.maxRange);
    lines = currentMap->lines;
  }
  
  vector2f locPrev = loc;
  float anglePrev = angle;
  float weight;
  assert (locCorrectionP0.size() > particleIdx);
  locCorrectionP0[particleIdx] = loc;
  for(int i=0; beingRefined && i<pointCloudParams.numSteps; i++){
    getPointCloudGradient(particleIdx, loc, angle, locGrad, angleGrad, pointCloud, pointNormals, weight, pointCloudParams, lineCorrespondences, lines);
    if(i==0) initialWeight = exp(weight);
    loc -= pointCloudParams.etaLoc*locGrad;
    angle -= pointCloudParams.etaAngle*angleGrad;
    beingRefined = fabs(angleGrad)>pointCloudParams.minRefineFraction*pointCloudParams.maxAngleGradient && locGrad.sqlength()>sq(pointCloudParams.minRefineFraction*pointCloudParams.maxLocGradient);
  }
  assert (locCorrectionP1.size() > particleIdx);
  locCorrectionP1[particleIdx] = loc;
  finalWeight = exp(weight);
  if(debug) printf("gradient: %7.3f,%7.3f %6.1f\u00b0\n",V2COMP(loc-locPrev),DEG(angle-anglePrev));
}

void VectorLocalization2D::refineLidar(const LidarParams &lidarParams)
{
  //FunctionTimer ft(__PRETTY_FUNCTION__);
  double tStart = GetTimeSec();
  laserEval.stage0Weights = 0.0;
  laserEval.stageRWeights = 0.0;
  laserEval.lastRunTime = GetTimeSec();
  
  // Transform laserpoints to robot frame
  vector< Vector2f > laserPoints(lidarParams.numRays);
  for(int i=0; i<lidarParams.numRays; i++){
    laserPoints[i] = lidarParams.laserToBaseTrans + lidarParams.laserToBaseRot*lidarParams.scanHeadings[i]*lidarParams.laserScan[i];
  }
  
  if(lidarParams.numSteps>0){
    for(int i=0; i<numParticles; i++){
      float initialWeight, finalWeight;
      refineLocationLidar(particles[i].loc, particles[i].angle, initialWeight, finalWeight, lidarParams, laserPoints);
      laserEval.stage0Weights += initialWeight;
      laserEval.stageRWeights += finalWeight;
    }
  }
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

void VectorLocalization2D::getParticles(std::vector<Particle2D> &_particles)
{
  _particles = particles;
}

void VectorLocalization2D::refinePointCloud(const vector<vector2f> &pointCloud, const vector<vector2f> &pointNormals, const PointCloudParams &pointCloudParams)
{
  //FunctionTimer ft(__PRETTY_FUNCTION__);
  double tStart = GetTimeSec();
  pointCloudEval.stage0Weights = 0.0;
  pointCloudEval.stageRWeights = 0.0;
  pointCloudEval.lastRunTime = GetTimeSec();
  
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
  localization->refineLocationPointCloud(v.id(), v.data().loc, v.data().angle,
      PARTICLE_POINT_CLOUD_EVAL[v.id()].stage0Weights, PARTICLE_POINT_CLOUD_EVAL[v.id()].stageRWeights,
      *REFINEMENT->pointCloud, *REFINEMENT->pointNormals, *REFINEMENT->pointCloudParams);

  particles[v.id()] = v.data();
}

void updatePointCloudParticle(graph_type::vertex_type& v)
{
  static const bool debug = true;

  if(debug) printf("\nParticle weights:\n");

  // Compute importance weights = observation x motion / samplingDensity
  float w1 = localization->observationWeightPointCloud(v.data().loc, v.data().angle, *UPDATE->pointCloud, *UPDATE->pointNormals, *UPDATE->pointCloudParams);
  float w2 = localization->motionModelWeight(v.data().loc, v.data().angle, *UPDATE->motionParams);
  v.data().weight = w1*w2/SAMPLING_DENSITY[v.id()]/TOTAL_DENSITY;

  if(debug) printf("%2lu: %f , %f , %f\n", v.id(), w1, w2, v.data().weight);
}

void samplingDensityPointCloudParticle(graph_type::vertex_type& v)
{
  static const bool debug = false;

  //Compute the sampling density
  float sqDensityKernelSize = sq(UPDATE->pointCloudParams->kernelSize);

  if(debug) printf("\nParticle samplingDensity:\n");

  float w = 0.99;
  for(int i=0; i<localization->getNumParticles(); i++){
    if(i==v.id())
      continue;

    if( (particles[i].loc - v.data().loc).sqlength() < sqDensityKernelSize && fabs(angle_diff(particles[i].angle, v.data().angle))<RAD(20.0))
      w++;
  }

  assert(SAMPLING_DENSITY.size() > v.id());
  SAMPLING_DENSITY[v.id()] = w;

  if(debug) printf("%2lu:%f\n", v.id(), w);
}

void VectorLocalization2D::computeLocation(vector2f& loc, float& angle)
{
  //sum all the poses (taking care od adding up headings appropriately)
  PoseReducer r = graph->map_reduce_vertices<PoseReducer>(PoseReducer::getPose);
  //maximum likelihood estimation of the pose
  currentLocation.x = r.x/numParticles;
  currentLocation.y = r.y/numParticles;
  currentAngle = std::atan2(r.heading_y, r.heading_x);

  currentLocStdDev.zero();
  currentAngleStdDev = 0.0;
  float xStdDev=0.0, yStdDev=0.0;
  for(int i=0; i<numParticles; i++){
    /*
    xStdDev += sq(float(particles[i].loc.x-currentLocation.x));
    yStdDev += sq(float(particles[i].loc.y-currentLocation.y));
    currentAngleStdDev += sq(float(angle_diff(particles[i].angle,currentAngle)));
    */
    xStdDev = max(xStdDev, float(fabs(particles[i].loc.x-currentLocation.x)));
    yStdDev = max(yStdDev, float(fabs(particles[i].loc.y-currentLocation.y)));
    currentAngleStdDev = max(currentAngleStdDev, float(fabs(angle_diff(particles[i].angle,currentAngle))));
  }
  /*
  currentAngleStdDev = sqrt(currentAngleStdDev/float(numParticles-1.0));
  xStdDev = sqrt(xStdDev/float(numParticles));
  yStdDev = sqrt(yStdDev/float(numParticles));
  */
  currentLocStdDev.set(xStdDev,yStdDev);
  loc = currentLocation;
  angle = currentAngle;
}

void VectorLocalization2D::resample(Resample type)
{
  switch(type){
    case NaiveResampling:{
      naiveResample();
    }break;
    case LowVarianceResampling:{
      lowVarianceResample();
    }break;
    case SensorResettingResampling:{
    }break;
  }
}

void VectorLocalization2D::lowVarianceResample()
{ 
  vector<Particle2D> newParticles;
  newParticles.resize(numParticles);
  float totalWeight = 0.0;
  float newWeight = 1.0/float(numParticles);
  int numRefinedParticles = (int) particles.size();
  
  refinedImportanceWeights = unrefinedImportanceWeights = 0.0;
  for(int i=0; i<numRefinedParticles; i++){
    //Get rid of particles with undefined weights
    if(isnan(particles[i].weight) || isinf(particles[i].weight) || particles[i].weight<0.0)
      particles[i].weight = 0.0;
    totalWeight += particles[i].weight;
    if(i<numParticles)
      refinedImportanceWeights += particles[i].weight;
    else
      unrefinedImportanceWeights += particles[i].weight;
  }
  
  if(totalWeight<FLT_MIN){
    //TerminalWarning("Particles have zero total weight!");
    for(int i=0; i<numParticles; i++){
      particles[i].weight = newWeight;
    }
    return;
    //exit(0);
  }
  
  float weightIncrement = totalWeight/float(numParticles);
  if(weightIncrement<FLT_MIN) TerminalWarning("Particle weights less than float precision");
  
  numRefinedParticlesSampled = numUnrefinedParticlesSampled = 0;
  float x = frand(0.0f,totalWeight);
  int j=0;
  float f=particles[0].weight;
  for(int i=0; i<numParticles; i++){
    while(f<x){
      j = (j+1)%numRefinedParticles;
      f += particles[j].weight;
    }
    if(j<numParticles)
      numRefinedParticlesSampled++;
    else
      numUnrefinedParticlesSampled++;
    
    newParticles[i] = particles[j];
    newParticles[i].weight = newWeight;
    if(particles[i].weight < FLT_MIN){
      //This particle was depleted: add replacement noise
      vector2f deltaLoc = vector2f(frand(-1.0,1.0),frand(-1.0,1.0))*0.05;
      float deltaAngle = frand(-1.0,1.0)*RAD(5.0);
      newParticles[i].loc += deltaLoc;
      newParticles[i].angle += deltaAngle;
    }
    x += weightIncrement;
  }
  
  particles = newParticles;
}

void VectorLocalization2D::naiveResample()
{
  vector<Particle2D> newParticles;
  static vector<Particle2D> oldParticles;
  newParticles.resize(numParticles);
  float totalWeight = 0.0;
  float newWeight = 1.0/numParticles;
  int numRefinedParticles = (int) particles.size();
  
  refinedImportanceWeights = unrefinedImportanceWeights = 0.0;
  int numInfs=0, numNans=0, numNegs=0;
  for(int i=0; i<numRefinedParticles; i++){
    if(isnan(particles[i].weight))
      numNans++;
    else if(isinf(particles[i].weight))
      numInfs++;
    else if(particles[i].weight<0.0)
      numNegs++;
    //Get rid of particles with undefined weights
    if(isnan(particles[i].weight) || isinf(particles[i].weight) || particles[i].weight<0.0)
      particles[i].weight = 0.0;
    totalWeight += particles[i].weight;
    if(i<numParticles)
      refinedImportanceWeights += particles[i].weight;
    else
      unrefinedImportanceWeights += particles[i].weight;
  }
  if(totalWeight<FLT_MIN){
    TerminalWarning("Particles have zero total weight!");
    printf("inf:%d nan:%d neg:%d\n",numInfs, numNans, numNegs);
    particles = oldParticles;
    return;
    //exit(0);
  }
  
  numRefinedParticlesSampled = numUnrefinedParticlesSampled = 0;
  for(int i=0; i<numParticles; i++){
    float x = frand(0.0f,totalWeight);
    float f=particles[0].weight;
    int j=0;
    while(f<x && j<numRefinedParticles-1){
      j++;
      f += particles[j].weight;
    }
    if(j<numParticles)
      numRefinedParticlesSampled++;
    else
      numUnrefinedParticlesSampled++;
    newParticles[i] = particles[j];
    newParticles[i].weight = newWeight;
  }
  particles = newParticles;
  oldParticles = particles;
}

void VectorLocalization2D::drawDisplay(vector<float> &lines_p1x, vector<float> &lines_p1y, vector<float> &lines_p2x, vector<float> &lines_p2y, vector<uint32_t> &lines_color,
                                     vector<float> &points_x, vector<float> &points_y, vector<uint32_t> &points_color, 
                                       vector<float> &circles_x, vector<float> &circles_y, vector<uint32_t> &circles_color, float scale)
{
  static const bool debug = false;
  static const bool drawParticles = true;
  static const bool drawLidarPoints = true;
  static const bool drawKinectPoints = true;
  static const bool drawCorrections = true;
  static const bool drawSceneRender = false;
  
  static const float GradientsScale = 1.0;
  static const int ParticleCloudColor = 0xF77274;
  static const int LidarPointColor = 0xF0761F;
  static const int LidarGradientColor = 0x2F95ED;
  static const int CloudPointColor = 0xDE2352;
  static const int CloudGradientColor = 0x8A23DE;
  static const int LineCorrespondenceColor = 0x93FA70;
  static const int LineExtentColor = 0xFFD659;
  static const int DebugColour = 0xFF0000;
  
  for(int i=0; i<debugLines.size(); i++){
    lines_p1x.push_back(scale*debugLines[i].P0().x);
    lines_p1y.push_back(scale*debugLines[i].P0().y);
    lines_p2x.push_back(scale*debugLines[i].P1().x);
    lines_p2y.push_back(scale*debugLines[i].P1().y);
    lines_color.push_back(DebugColour);
  }
  debugLines.clear();
  
  if(drawCorrections){
    assert(locCorrectionP0.size() == numParticles);
    assert(locCorrectionP1.size() == numParticles);

    for(int i=0; i<locCorrectionP0.size(); i++){
      vector2f p2 = locCorrectionP1[i] + 20.0*(locCorrectionP1[i] - locCorrectionP0[i]);
      lines_p1x.push_back(scale*locCorrectionP0[i].x);
      lines_p1y.push_back(scale*locCorrectionP0[i].y);
      lines_p2x.push_back(scale*p2.x);
      lines_p2y.push_back(scale*p2.y);
      lines_color.push_back(0x1BE042);
    }
  }
  
  if(drawParticles){
    //Draw all (unrefined) particles
    for(int i=0; i<particles.size(); i++){
      circles_x.push_back(scale*particles[i].loc.x);
      circles_y.push_back(scale*particles[i].loc.y);
      circles_color.push_back(ParticleCloudColor);
      
      vector2f h;
      h.heading(particles[i].angle);
      vector2f p = particles[i].loc + 0.3*h;
      lines_p1x.push_back(scale*particles[i].loc.x);
      lines_p1y.push_back(scale*particles[i].loc.y);
      lines_p2x.push_back(scale*p.x);
      lines_p2y.push_back(scale*p.y);
      lines_color.push_back(ParticleCloudColor);
    }
  }
  
  if(drawSceneRender){
    vector2f loc;
    float angle;
    computeLocation(loc,angle);
    vector<line2f> lines = currentMap->sceneRender(loc,0.0,RAD(359.9));
    vector<int> locVisibilityList;
    if(currentMap->preRenderExists)
      locVisibilityList = *currentMap->getVisibilityList(loc.x, loc.y);
    else
      locVisibilityList = currentMap->getSceneLines(loc, 6.5);
    
    for(unsigned int i=0; i<locVisibilityList.size(); i++){
      lines_p1x.push_back(scale*currentMap->lines[locVisibilityList[i]].P0().x);
      lines_p1y.push_back(scale*currentMap->lines[locVisibilityList[i]].P0().y);
      lines_p2x.push_back(scale*currentMap->lines[locVisibilityList[i]].P1().x);
      lines_p2y.push_back(scale*currentMap->lines[locVisibilityList[i]].P1().y);
      lines_color.push_back(0x00DE07);
    }
    
    vector<vector2f> rays;
    vector2f ray;
    bool duplicate;
    for(unsigned int i=0; i<lines.size(); i++){
      lines_p1x.push_back(scale*lines[i].P0().x);
      lines_p1y.push_back(scale*lines[i].P0().y);
      lines_p2x.push_back(scale*lines[i].P1().x);
      lines_p2y.push_back(scale*lines[i].P1().y);
      lines_color.push_back(DebugColour);
      
      duplicate = false;
      ray = lines[i].P0()-loc;
      for(int j=0; j<rays.size()&&!duplicate; j++){
        if(ray.norm().dot(rays[j].norm())>cos(RAD(0.1)) && fabs(ray.length()-rays[j].length())<0.01){
          duplicate = true;
          j = rays.erase(rays.begin()+j) - rays.begin();
        }
      }
      if(!duplicate)
        rays.push_back(ray);
      
      duplicate = false;
      ray = lines[i].P1()-loc;
      for(int j=0; j<rays.size()&&!duplicate; j++){
        if(ray.norm().dot(rays[j].norm())>cos(RAD(0.1)) && fabs(ray.length()-rays[j].length())<0.01){
          duplicate = true;
          j = rays.erase(rays.begin()+j) - rays.begin();
        }
      }
      if(!duplicate)
        rays.push_back(ray);
    }
    
    for(int i=0; i<rays.size(); i++){
      vector2f p0 = loc;
      vector2f p1 = rays[i]+loc;
      lines_p1x.push_back(scale*p0.x);
      lines_p1y.push_back(scale*p0.y);
      lines_p2x.push_back(scale*p1.x);
      lines_p2y.push_back(scale*p1.y);
      lines_color.push_back(LineExtentColor);
    }
    
  }
  
  if(drawLidarPoints){
    // Draw all LIDAR points considered for gradient descent, as well as their gradients
    for(int i=0; i<points.size(); i++){
      points_x.push_back(scale*points[i].x());
      points_y.push_back(scale*points[i].y());
      points_color.push_back(LidarPointColor);
      Vector2f p = points[i] - GradientsScale*gradients[i];
      lines_p1x.push_back(scale*points[i].x());
      lines_p1y.push_back(scale*points[i].y());
      lines_p2x.push_back(scale*p.x());
      lines_p2y.push_back(scale*p.y());
      lines_color.push_back(LidarGradientColor);
    }
  }
  
  if(drawKinectPoints){
    assert(POINTS2.size() == GRADIENTS2.size());
    assert(POINTS2.size() == numParticles);

    // Draw all point cloud points considered for gradient descent, as well as their gradients
    for (int particleIdx=0; particleIdx<POINTS2.size(); particleIdx++){
      const std::vector<Vector2f>& points2 = POINTS2[particleIdx];
      const std::vector<Vector2f>& gradients2 = GRADIENTS2[particleIdx];
      assert(points2.size() == gradients2.size());

      for(int i=0; i<points2.size(); i++){
        points_x.push_back(scale*points2[i].x());
        points_y.push_back(scale*points2[i].y());
        points_color.push_back(CloudPointColor);
        Vector2f p = points2[i] - GradientsScale*gradients2[i];
        lines_p1x.push_back(scale*points2[i].x());
        lines_p1y.push_back(scale*points2[i].y());
        lines_p2x.push_back(scale*p.x());
        lines_p2y.push_back(scale*p.y());
        lines_color.push_back(CloudGradientColor);
      }
    }
  }
  
}

void VectorLocalization2D::saveProfilingStats(FILE* f)
{
  fprintf(f, "%f, %f, ",refineTime, updateTime);
}

void VectorLocalization2D::saveRunLog(FILE* f)
{
  static const bool saveParticles = false;
  fprintf(f, "%d, %d, %.15f, %.15f, ",numRefinedParticlesSampled, numUnrefinedParticlesSampled, refinedImportanceWeights, unrefinedImportanceWeights);
  if(saveParticles){
    for(int i=0; i<particles.size(); i++){
      fprintf(f, "%.3f, %.3f, %.2f, %f, ",V2COMP(particles[i].loc), DEG(particles[i].angle), particles[i].weight);
    }
  }
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

SamplingDensityReducer SamplingDensityReducer::computeTotalDensity(const graph_type::vertex_type& v)
{
  SamplingDensityReducer r;
  r.samplingDensity = SAMPLING_DENSITY[v.id()];
  return r;
}

SamplingDensityReducer& SamplingDensityReducer::operator+=(const SamplingDensityReducer& other)
{
  samplingDensity += other.samplingDensity;
  return *this;
}

void initializeParticle(graph_type::vertex_type& v)
{
  v.data() = Particle2D(
      randn(PARTICLE_INITIALIZER->locationUncertainty, PARTICLE_INITIALIZER->loc.x),
      randn(PARTICLE_INITIALIZER->locationUncertainty, PARTICLE_INITIALIZER->loc.y),
      randn(PARTICLE_INITIALIZER->angleUncertainty, PARTICLE_INITIALIZER->angle),
      1.0);

  particles[v.id()] = v.data();
}

VectorLocalization2D::ParticleInitializer::ParticleInitializer()
{
  loc.zero();
  angle = 0.0;
  locationUncertainty = 0.0;
  angleUncertainty = 0.0;
}

VectorLocalization2D::ParticleInitializer::ParticleInitializer(vector2f _loc, float _angle, float _locationUncertainty, float _angleUncertainty)
{
  loc.x = _loc.x;
  loc.y = _loc.y;
  angle = _angle;
  locationUncertainty = _locationUncertainty;
  angleUncertainty = _angleUncertainty;
}

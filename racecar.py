import math
import random
import pylab

class Coordinate(object):
    """
    A Coordinate represents a point in the cartesian x,y,z space.
    """
    
    def __init__(self,x,y,z=0):
        """
        Initialize with x, y, and z points. Z initializes to 0 by default.
        """
        self.x = x
        self.y = y
        self.z = z
    
    def getX(self):
        return self.x
        
    def getY(self):
        return self.y
        
    def getZ(self):
        return self.z
        
    def __str__(self):
        return "(%0.2f, %0.2f, %0.2f)" % (self.x, self.y, self.z)

class Turn(object):
    """
    A Turn is an arc representing the centreline of a single turn of a Track.
    centre:    centre point of the circle
    radius:    radius of the circle
    theta1:    start of arc
    theta2:    end of arc
    length:    length of arc
    """
    
    def __init__(self, centre, radius, thetaStart=0, thetaEnd=2*math.pi):
        self.centre = centre
        self.radius = radius
        self.thetaStart = thetaStart
        self.thetaEnd = thetaEnd
        
    def getCentre(self):
        return self.centre
        
    def getRadius(self):
        return self.radius
        
    def getThetaStart(self):
        return self.thetaStart
        
    def getThetaEnd(self):
        return self.thetaEnd
        
    def getLength(self):
        """
        Calculates the length of the turn arc.
        """
        if self.thetaStart < self.thetaEnd:
            self.length = 2*math.pi*self.radius*(abs(self.thetaStart-self.thetaEnd))/(2*math.pi)
        else:
            self.length = 2*math.pi*self.radius*(2*math.pi-abs(self.thetaStart-self.thetaEnd))/(2*math.pi)
        return self.length
        
    def getPoints(self,direction="left",distanceBetweenPoints=10):
        """
        Gets points from thetaStart to thetaEnd, counterclockwise if turn direction is left, clockwise if right.
        """
        self.points = []
        if direction == "right":
            start = self.thetaStart
            end = self.thetaEnd
            self.thetaStart = end
            self.thetaEnd = start
        if self.thetaStart < self.thetaEnd:
            deltaTheta = abs(self.thetaStart-self.thetaEnd)/(self.getLength()/distanceBetweenPoints)
        else:
            deltaTheta = (2*math.pi-abs(self.thetaStart-self.thetaEnd))/(self.getLength()/distanceBetweenPoints)
        for i in range(int(self.getLength()/distanceBetweenPoints)):
            theta = self.thetaStart+i*deltaTheta
            x = self.centre.getX()+self.getRadius()*math.cos(theta)
            y = self.centre.getY()+self.getRadius()*math.sin(theta)
            z = self.centre.getZ()+0
            self.points.append(Coordinate(x,y,z))
        return self.points
        
    def __str__(self):
        return "("+str(self.centre)+", r=%0.2f, start angle=%0.2f, end angle=%0.2f, length=%0.2f)" % (self.radius,self.thetaStart,self.thetaEnd,self.getLength())

class Straight(object):
    """
    A Straight is an line representing the centreline of a straight between two
    Turns on a Track.
    start:    starting point of the straight
    end:      ending point of the straight
    length:   length of the straight
    theta:    angle to line
    """
    
    def __init__(self, start, end):
        self.start = start
        self.end = end
        
    def getStart(self):
        return self.start
        
    def getEnd(self):
        return self.end
        
    def getLength(self):
        xStart = self.start.getX()
        yStart = self.start.getY()
        xEnd = self.end.getX()
        yEnd = self.end.getY()
        self.length = math.sqrt((xStart-xEnd)**2+(yStart-yEnd)**2)
        return self.length
    
    def getTheta(self):
        adjacent = self.end.getX()-self.start.getX()
        hypotenuse = self.getLength()
        if self.end.getY() >= self.start.getY():
            self.theta = math.acos(adjacent / self.getLength())
        else:
            self.theta = -math.acos(adjacent / self.getLength())
        return self.theta
        
    def getPoints(self,distanceBetweenPoints=10):
        self.points = []
        deltaX = distanceBetweenPoints * math.cos(self.getTheta())
        deltaY = distanceBetweenPoints * math.sin(self.getTheta())
        for i in range(int(self.getLength()/distanceBetweenPoints)):
            x = self.start.getX()+i*deltaX
            y = self.start.getY()+i*deltaY
            z = self.start.getZ()+0
            self.points.append(Coordinate(x,y,z))
        return self.points
    
    def onSegment(self,p,q,r):
        if(q.x <= max(p.x, r.x) and q.x >= min(p.x, r.x) and q.y <= max(p.y, r.y) and q.y >= min(p.y, r.y)):
            return True
        return False
        
    def orientation(self,p,q,r):
        val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
        if (val == 0):
            return 0
        elif (val > 0):
            return 1
        else:
            return 2
    
    def intersect(self, straight):
        p1 = self.getStart()
        q1 = self.getEnd()
        p2 = straight.getStart()
        q2 = straight.getEnd()
        o1 = self.orientation(p1, q1, p2)
        o2 = self.orientation(p1, q1, q2)
        o3 = self.orientation(p2, q2, p1)
        o4 = self.orientation(p2, q2, q1)
        if (o1 != o2 and o3 != o4):
            return True
 
    # Special Cases
    # p1, q1 and p2 are colinear and p2 lies on segment p1q1
        if (o1 == 0 and self.onSegment(p1, p2, q1)):
            return True
 
    # p1, q1 and p2 are colinear and q2 lies on segment p1q1
        if (o2 == 0 and self.onSegment(p1, q2, q1)):
            return True
 
    # p2, q2 and p1 are colinear and p1 lies on segment p2q2
        if (o3 == 0 and self.onSegment(p2, p1, q2)):
            return True
 
    # p2, q2 and q1 are colinear and q1 lies on segment p2q2
        if (o4 == 0 and self.onSegment(p2, q1, q2)):
            return True
            
    # none of the above cases apply
        return False        
        
    def __str__(self):
        return "("+str(self.start)+", angle=%0.2f, length=%0.2f)" % (self.getTheta(),self.getLength())

class SimpleTrack(object):
    """
    A Simple Track represents the centreline of a racetrack.
    """
    
    def __init__(self,trackRadius,numTurns):
        """
        Initializes a Simple Track using a number of straights connected by
        arcs. Input is overall track radius (metres) and number of turns.
        
        turns:      list of turns that define the track
        theta:      rotational angle from centre of track to centre of turn
        distance:   radial distance from centre of track to centre of turn
        
        trackCircumferencePoints:   conceptual construction of track circumference
        """
        
        self.trackCentre = Coordinate(0.0,0.0)
        self.trackRadius = trackRadius
        self.numTurns = numTurns
        self.trackCircumferencePoints = []
        self.turns = []
        self.straights = []
        self.curves = []
        self.curveDirections = []
    
    def getRadius(self):
        return self.trackRadius
        
    def getNumTurns(self):
        return self.numTurns
    
    def getCircumference(self):        
        """
        Establishes track circumference.
        """
        self.trackCircumference = Turn(self.trackCentre, self.getRadius())
        return self.trackCircumference
    
    def getTurns(self):
        """
        Establishes track turns sequentially around the track circumference and
        saves them to a list of turns. For each turn, establishes parameters,
        and then checks to make sure newly created turn does not overlap any
        existing ones.
        """
        
        for num in range(self.numTurns):
            turns_overlap = True
            while turns_overlap == True:
                #Establish turn parameters 
                theta = 2*math.pi * (float(num)+random.random()) / float(self.numTurns)
                distance = random.gauss(self.trackRadius, float(self.trackRadius)/5.0)
                x = distance * math.cos(theta)
                y = distance * math.sin(theta)
                z = 0
                c = Coordinate(x,y,z)
                r = (5.0/float(self.numTurns))*random.uniform(0.05,1.0)*float(self.trackRadius)
                turn = Turn(c,r)
                            
                #Avoid turn circumference overlap
                distances_between_turns = []
                if num == 0:
                    turns_overlap = False
                else:
                    for turn in self.turns:
                        xi = turn.getCentre().getX()
                        yi = turn.getCentre().getY()
                        ri = turn.getRadius()
                        distance_between_centres = math.sqrt((xi-x)**2+(yi-y)**2)
                        distance_between_turns = distance_between_centres - (ri+r)
                        distances_between_turns.append(distance_between_turns)
                    if any(distance < 0 for distance in distances_between_turns):
                        turns_overlap = True
                    else:
                        turns_overlap = False
                    
            #Adds turn to turns list
            self.turns.append(Turn(c, r))
        return self.turns
        
    def getTangents(self):
        """
        Add Tangents between Turns. Identifies both external tangents and both
        internal tangents between each sequential pair of Turns.
        """
        
        self.tangents = []
        
        if self.turns == []:
            self.getTurns()
        
        for i in range(self.numTurns):
            if self.numTurns == 2 and i == 2:
                break
            tangent = []
            turn1Radius = self.turns[i-1].getRadius()
            turn2Radius = self.turns[i].getRadius()
            if turn1Radius > turn2Radius:
                c1 = self.turns[i].getCentre()
                c2 = self.turns[i-1].getCentre()
                r1 = self.turns[i].getRadius()
                r2 = self.turns[i-1].getRadius()
            else:
                c1 = self.turns[i-1].getCentre()
                c2 = self.turns[i].getCentre()
                r1 = turn1Radius
                r2 = turn2Radius
            c1x = c1.getX()
            c1y = c1.getY()
            c2x = c2.getX()
            c2y = c2.getY()
            
            #External tangents
            D = math.sqrt((c1x-c2x)**2+(c1y-c2y)**2)
            d = r1*D/(r2-r1)
            l = math.sqrt(d**2-r1**2)
            L = math.sqrt((D+d)**2-r2**2)
            adjacent = c1x-c2x
            hypotenuse = D
            theta = math.acos(adjacent / hypotenuse)
            alpha = math.asin(r1/d)
            beta = theta - alpha
            gamma = theta + alpha
            if c1y >= c2y:
                PE = Coordinate(c1x+d*math.cos(theta),c1y+d*math.sin(theta))
                T1Ea = Coordinate(PE.getX()-l*math.cos(beta),PE.getY()-l*math.sin(beta))
                T2Ea = Coordinate(PE.getX()-L*math.cos(beta),PE.getY()-L*math.sin(beta))
                T1Eb = Coordinate(PE.getX()-l*math.cos(gamma),PE.getY()-l*math.sin(gamma))
                T2Eb = Coordinate(PE.getX()-L*math.cos(gamma),PE.getY()-L*math.sin(gamma))
            else:
                PE = Coordinate(c1x+d*math.cos(theta),c1y-d*math.sin(theta))
                T1Ea = Coordinate(PE.getX()-l*math.cos(beta),PE.getY()+l*math.sin(beta))
                T2Ea = Coordinate(PE.getX()-L*math.cos(beta),PE.getY()+L*math.sin(beta))
                T1Eb = Coordinate(PE.getX()-l*math.cos(gamma),PE.getY()+l*math.sin(gamma))
                T2Eb = Coordinate(PE.getX()-L*math.cos(gamma),PE.getY()+L*math.sin(gamma))
            if turn1Radius > turn2Radius:
                TEa = Straight(T2Ea,T1Ea)
                TEb = Straight(T2Eb,T1Eb)
            else:
                TEa = Straight(T1Ea,T2Ea)
                TEb = Straight(T1Eb,T2Eb)

            #Internal tangents
            #D is the same
            d = r1*D/(r2+r1)
            l = math.sqrt(d**2-r1**2)
            L = math.sqrt((D-d)**2-r2**2)
            #theta is the same
            alpha = math.asin(r1/d)
            beta = theta - alpha
            gamma = theta + alpha
            if c1y >= c2y:
                PI = Coordinate(c1x-d*math.cos(theta),c1y-d*math.sin(theta))
                T1Ia = Coordinate(PI.getX()+l*math.cos(gamma),PI.getY()+l*math.sin(gamma))
                T2Ia = Coordinate(PI.getX()-L*math.cos(gamma),PI.getY()-L*math.sin(gamma))
                T1Ib = Coordinate(PI.getX()+l*math.cos(beta),PI.getY()+l*math.sin(beta))
                T2Ib = Coordinate(PI.getX()-L*math.cos(beta),PI.getY()-L*math.sin(beta))
            else:
                PI = Coordinate(c1x-d*math.cos(theta),c1y+d*math.sin(theta))
                T1Ia = Coordinate(PI.getX()+l*math.cos(gamma),PI.getY()-l*math.sin(gamma))
                T2Ia = Coordinate(PI.getX()-L*math.cos(gamma),PI.getY()+L*math.sin(gamma))
                T1Ib = Coordinate(PI.getX()+l*math.cos(beta),PI.getY()-l*math.sin(beta))
                T2Ib = Coordinate(PI.getX()-L*math.cos(beta),PI.getY()+L*math.sin(beta))
            if turn1Radius > turn2Radius:
                TIa = Straight(T2Ia,T1Ia)
                TIb = Straight(T2Ib,T1Ib)
            else:
                TIa = Straight(T1Ia,T2Ia)
                TIb = Straight(T1Ib,T2Ib)
            
            #add tangents to tangents list
            tangent.append(TEa)
            tangent.append(TEb)
            tangent.append(TIa)
            tangent.append(TIb)
            self.tangents.append(tangent)
            
        return self.tangents
    
    def getStraights(self):
        for i in range(len(self.tangents)):
            tangentSet = self.tangents[i]
            TEaStart = tangentSet[0].getStart()
            TEbStart = tangentSet[1].getStart()
            TEaDist = Straight(self.trackCentre,TEaStart).getLength()
            TEbDist = Straight(self.trackCentre,TEbStart).getLength()
            if TEaDist > TEbDist:
                straight = tangentSet[0]
            else:
                straight = tangentSet[1]
            self.straights.append(straight)
            self.curveDirections.append(str("left"))
            
        for i in range(len(self.straights)):
            tangentSet1 = self.tangents[i-1]
            tangentSet2 = self.tangents[i]
            
            #check if straights i-1 and i intersect
            if self.straights[i-1].intersect(self.straights[i]):
                self.curveDirections[i] = str("right")
                TIaStart = tangentSet1[2].getStart()
                TIbStart = tangentSet1[3].getStart()
                TIaDist = Straight(self.trackCentre,TIaStart).getLength()
                TIbDist = Straight(self.trackCentre,TIbStart).getLength()
                if TIaDist > TIbDist:
                    self.straights[i-1] = self.tangents[i-1][2]
                else:
                    self.straights[i-1] = self.tangents[i-1][3]
                    
                TIaStart = tangentSet2[2].getStart()
                TIbStart = tangentSet2[3].getStart()
                TIaDist = Straight(self.trackCentre,TIaStart).getLength()
                TIbDist = Straight(self.trackCentre,TIbStart).getLength()
                if TIaDist > TIbDist:
                    self.straights[i] = self.tangents[i][3]
                else:
                    self.straights[i] = self.tangents[i][2]    
            
        return self.straights
    
    def getCurves(self):
        for i in range(len(self.turns)):
            centre = self.turns[i-1].getCentre()
            radius = self.turns[i-1].getRadius()
            startPoint = self.straights[i-1].getEnd()
            endPoint = self.straights[i].getStart()
            
            if startPoint.x > centre.x:
                thetaStart = math.asin((startPoint.y-centre.y)/radius)
            else:
                thetaStart = math.pi - math.asin((startPoint.y-centre.y)/radius)
                
            if endPoint.x > centre.x:
                thetaEnd = math.asin((endPoint.y-centre.y)/radius)
            else:
                thetaEnd = math.pi - math.asin((endPoint.y-centre.y)/radius)
                
            curve = Turn(centre,radius,thetaStart,thetaEnd)
            self.curves.append(curve)
        return self.curves
    
    def plotEverything(self):
        """
        #Plot everything
        """
        font = {'family' : 'normal', 'size' : 24}
        pylab.rc('font',**font)
        circumference = self.getCircumference()
        points = circumference.getPoints()
        centre = circumference.getCentre()
        #pylab.plot(centre.getX(),centre.getY(),'bo')
        #for point in points:
        #    pylab.plot(point.getX(),point.getY(),'bo')
        """
        for turn in self.turns:
            points = turn.getPoints()
            centre = turn.getCentre()
            for point in points:
                pylab.plot(point.getX(),point.getY(),'ro')
            pylab.plot(centre.getX(),centre.getY(),'ro')
            index = self.turns.index(turn)
            #if index == 1:
            #    pylab.plot(centre.getX(),centre.getY(),'go',markersize=10)
            #    break
        """
        for straight in self.straights:
            points = straight.getPoints()
            for point in points:
                pylab.plot(point.getX(),point.getY(),'bo')
            startPoint = straight.getStart()
            pylab.plot(startPoint.getX(),startPoint.getY(),'ro',markersize=20)
            endPoint = straight.getEnd()
            pylab.plot(endPoint.getX(),endPoint.getY(),'go',markersize=20)
            index = self.straights.index(straight)
            #if index == 1:
            #    break
                
        for curve in self.curves:
            i = self.curves.index(curve)
            print curve
            points = curve.getPoints(self.curveDirections[i])
            for point in points:
                pylab.plot(point.getX(),point.getY(),'go',markersize=10)
        """
        for tangentSet in self.tangents:
            for tangent in tangentSet:
                points = tangent.getPoints()
                for point in points:
                    pylab.plot(point.getX(),point.getY(),'bo')
        """
        """
        track_min_dimension = 0
        for i in range(len(turns_circumferences_x)):
            pylab.plot(turns_circumferences_x[i],turns_circumferences_y[i],'bo')
            min_dimension = min(min(turns_circumferences_x[i]),min(turns_circumferences_y[i]),-(max(turns_circumferences_x[i])),-(max(turns_circumferences_y[i]))) 
            if min_dimension < track_min_dimension:
                track_min_dimension = min_dimension
        for i in range(len(turns)):
            pylab.plot(P1[i][0],P1[i][1],'ro',markersize=20)
        for i in range(len(turns_centre_x)):
            pylab.annotate(str(i),(turns_centre_x[i],turns_centre_y[i]))    
        #pylab.xlim(track_min_dimension,-track_min_dimension)
        #pylab.ylim(track_min_dimension,-track_min_dimension)
        """
def runSim(radius,numTurns):
    track = SimpleTrack(radius,numTurns)
    track.getTurns()
    track.getTangents()
    track.getStraights()
    track.getCurves()
    track.plotEverything()
    print track.curveDirections
    
def makeArc(centre = Coordinate(0.0,0.0), radius=100, thetaStart=0, thetaEnd=math.pi*2):
    turn2 = Turn(centre,radius)
    points2 = turn2.getPoints()
    for point in points2:
        pylab.plot(point.getX(),point.getY(),'ro')
    
    turn = Turn(centre,radius,thetaStart,thetaEnd)
    points = turn.getPoints()
    for point in points:
        pylab.plot(point.getX(),point.getY(),'go',markersize=10)
    pylab.plot(points[0].getX(),points[0].getY(),'bo',markersize=20)
    pylab.plot(points[-1].getX(),points[-1].getY(),'ro',markersize=20)
    print turn
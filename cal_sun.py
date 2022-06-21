import numpy as N

class SunPosition:

    def __init__(self):
        pass


    def days(self, dd, mm):
        '''
        Arguement:
        dd - int, the day in the month
        mm - str, the month

        reference: J Duffie pp14, Table 1.6.1
        '''

        if mm=='Jan':
            days=dd

        elif mm=='Feb':
            days=31+dd

        elif mm=='Mar':
            days=59+dd

        if mm=='Apr':
            days=90+dd

        elif mm=='May':
            days=120+dd

        elif mm=='Jun':
            days=151+dd

        if mm=='Jul':
            days=181+dd

        elif mm=='Aug':
            days=212+dd

        elif mm=='Sep':
            days=243+dd

        if mm=='Oct':
            days=273+dd

        elif mm=='Nov':
            days=304+dd

        elif mm=='Dec':
            days=334+dd

        return days


    def declination(self, days, form='detail'):
        '''
        Reference: Solar Engineering of Thermal Processes, 4th edition, John A. Duffie and William A. Beckman, page 13

        declination angle: delta=23.45*sin(360*(284+day)/365)
 
        Arguement:
        day - int, day of the year (1-365)
        form - str, 'detail' or simple' model

        Return:
        delta - declination angle (deg)
        '''
        if form=='detail':
            B=float(days-1)*360./365.*N.pi/180.

            delta=(180./N.pi)*(0.006918 - 0.399912*N.cos(B) +0.070257*N.sin(B)- 0.006758*N.cos(2.*B) + 0.000907*N.sin(2.*B)- 0.002697*N.cos(3.*B) + 0.00148*N.sin(3.*B))

        else:
            delta=23.45*N.sin(360.*float(284+days)/365.*N.pi/180.) # deg

        return delta


    def solarhour(self, delta, latitude):

        '''        
        Reference: Solar Engineering of Thermal Processes, 4th edition, John A. Duffie and William A. Beckman, page 17

        Arguement:
        delta: declination angle- float, deg
        latitude: latitude angle: float, deg

        return: 
        hour: length of the daylight hour
        sunrise: the solar hour angle of the sunrise, deg

        '''

        sunset=N.arccos(-N.tan(latitude*N.pi/180.)*N.tan(delta*N.pi/180.))*180./N.pi # deg

        sunrise=-sunset

        hour=(sunset-sunrise)/15.
       
        return hour, sunrise


    def zenith(self, latitude, delta, omega):
        '''
        ref. eq.1.6.5
        Arguement:
        latitude: latitude angle, float, deg
        delta:  declination angle, float, deg
        omega: solar hour angle, float, deg

        return:
        theta: the zenith angle, float, deg

        '''        
        
        latitude*=N.pi/180.
        delta*=N.pi/180.
        omega*=N.pi/180.

        theta=N.arccos(N.cos(latitude)*N.cos(delta)*N.cos(omega)+N.sin(latitude)*N.sin(delta))*180./N.pi

        return theta
        
    def azimuth(self, latitude, theta, delta, omega):

        '''
        ref: eq. 1.6.6
        from South to West
        Arguement:
        latitude: latitude angle, deg
        delta: declination angle ,deg
        theta: zenith angle, deg
        omega: solar hour angle, deg

        return:
        phi: azimuth angle, deg, from South to west
        '''
        latitude*=N.pi/180.
        delta*=N.pi/180.
        theta*=N.pi/180.
    
        a1=N.cos(theta)*N.sin(latitude)-N.sin(delta)
        a2=N.sin(theta)*N.cos(latitude)
        b=a1/a2

        if abs(b+1.)<1e-10:
            phi=N.pi

        elif abs(b-1.)<1e-10:
            phi=0.
        else:
            phi=abs(N.arccos((N.cos(theta)*N.sin(latitude)-N.sin(delta))/(N.sin(theta)*N.cos(latitude)))) # unit radian

        if omega<0:
            phi=-phi

        phi*=180./N.pi

        return phi

    def convert_AZEL_to_declination_hour(self, theta, phi, latitude):

        '''
        Arguement:
        theta: zenith angle, deg
        phi: azimuth angle deg
        latitude: latitude latitude , deg

        return:
        delta: declination angle, deg
        omega: solar hour angle, deg
        '''      
        phi*=N.pi/180.
        theta*=N.pi/180.
        latitude*=N.pi/180.

        delta=N.arcsin(N.cos(theta)*N.sin(latitude)-N.cos(abs(phi))*N.sin(theta)*N.cos(latitude))

        omega=N.arccos((N.cos(theta)-N.sin(latitude)*N.sin(delta))/(N.cos(latitude)*N.cos(delta)))
        if phi<0:
            omega=-omega

        delta*=180./N.pi
        omega*=180./N.pi

        return delta, omega



if __name__=='__main__':
	# example: PS10, summer solstice, solar noon
	latitude=34.96
	
	sun=SunPosition()
	dd=sun.days(22, 'Jun')
	delta=sun.declination(dd)
	print 'Declination angle', delta


	daytime,sunrise=sun.solarhour(delta, latitude)
	print 'Day timeS', daytime
	print 'sun rise', sunrise
	print ''
	
	
	omega=3.*15 # solar noon
	theta=sun.zenith(latitude, delta, omega)
	phi=sun.azimuth(latitude, theta, delta, omega)
	print 'elevation', 90.-theta
	print 'azimuth', phi

	



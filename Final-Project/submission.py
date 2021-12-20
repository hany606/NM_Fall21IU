# Microscopic model of traffic flow
# Sources:
# - https://gafferongames.com/post/integration_basics/
# - https://towardsdatascience.com/simulating-traffic-flow-in-python-ee1eab4dd20f
# Run it: python3 submission.py

from abc import abstractmethod
from math import sqrt, isclose
from time import sleep, time

class Physics:
    def __init__(self, x0=0.0, v0=0.0, a0=0.0):
        self.x = x0    # starting point
        self.v = v0
        self.a = a0


class Vehicle:
    def __init__(self, x0=0.0, v0=0.0, a0=0.0, id=None):
        self.p = Physics(x0, v0, a0)
        self.id = id

    def _compute_c_coeff(self, dxm):
        c = 0.0
        if(dxm >= 45):  # Free-flow
            c = 1075.0/dxm
        elif(dxm >= 25):   # Synchronized-flow
            c = 0.6*dxm + 362.0/dxm
        elif(dxm >= 7): # Jam-flow
            c = 0.6*dxm
        return c

    @abstractmethod
    def compute_dx(pm, pm1, len_traffic, track_length):
        # compute the delta x (difference between the cars)
        dx = pm1.x - pm.x
        # if it is single vehicles in the track
        # Or if the delta x is negative, in the case of the leader car in the starting of the road and the follower at the end
        #   Thus, we add the track_length to the delta x as it is ring road
        #   Kinda the same as taking the small angle instead of the big angle between two lines
        #   Example: dx should be 30 in case of (track_length=900): 10, 880 (dx=-870) -> 10+900, 880 (dx=30)
        if(dx < 0 or len_traffic == 1):
            dx += track_length
        return dx
    
    def update_acc(self, pm, pm1, len_traffic, track_length):
        self.dx = Vehicle.compute_dx(pm, pm1, len_traffic, track_length)
        # if dx is equal to zero, we will have division by zero, thus we avoid such a situation
        if not isclose(self.dx, 0, abs_tol=1e-10):
            c = self._compute_c_coeff(self.dx)
            dv = pm1.v - pm.v
            self.p.a = c*dv/self.dx
        self.p.a = 0 if isclose(abs(self.p.a), 0, abs_tol=1e-10) else self.p.a

    def update(self, dt, eps, track_length):
        # Taylor-series expansion
        # "Since this is only an approximation, the speed can become negative at times (but the model does not allow for that). An instability arises when the speed is negative, and the position and speed diverge into negative infinity"
        # Source: https://towardsdatascience.com/simulating-traffic-flow-in-python-ee1eab4dd20f
        # "To overcome this problem, whenever we predict a negative speed we will set it equal to zero and work out way from there:"
        if(self.p.v + self.p.a*dt < 0):
            self.p.x -= 0.5*self.p.v*self.p.v/self.p.a
            self.p.v = 0.0
        else:
            self.p.x += self.p.v*dt + self.p.a*dt*dt/2
            self.p.v += self.p.a*dt
        
        # Ring round track -> limit the position
        if(self.p.x >= track_length):
            self.p.x -= (self.p.x // track_length) * track_length

# Discrete microscopic traffic modeling
class Traffic:
    def __init__(self, track_length=None, tau=5, initial_velocity=40.0, initial_acceleration=0.0, epsilon=0.0001):
        self._set_parameters(track_length=track_length, tau=tau, initial_velocity=initial_velocity, initial_acceleration=initial_acceleration, epsilon=epsilon)
        self.reset()

    def _set_parameters(self, *args, **kwargs):
        self.track_length = kwargs["track_length"]
        self.initial_velocity = kwargs["initial_velocity"]
        self.initial_acceleration = kwargs["initial_acceleration"]
        self.tau = kwargs["tau"]    # time interval
        self.epsilon = kwargs["epsilon"]    # relative speed calculation error

    def reset(self):
        self.id = int(9e9)
        self.dt = self.tau
        # At the initial moment of time: 𝑡 = 0 [𝑠]
        #     the road is empty then at the point of entry: 𝑥 = 0 [m] one vehicle begins to appear every time
        #     𝜏 > 5 [𝑠] moving with initial speed 𝑣𝑚 = 40 [m/s] and zero acceleration 𝑣𝑚′ = 0.
        self.vehicles = [Vehicle(x0=0.0, v0=self.initial_velocity, a0=self.initial_acceleration, id=self.id)]
        self.timesteps = 0.0
        self.add_flag = False
        self.num_steps = 0

    def _compute_termination(self):
        # Method1: Working good  (The maximum of the velocity is less than or equal to the epsilon means it is termination) 
        # mx_val = -999e999
        # for i in range(len(self.vehicles)):
        #     mx_val = max(mx_val, self.vehicles[i].p.v)
        # # print(mx_val)
        # if(mx_val > self.epsilon):
        #     return False
        # else:
        #     return True
        # Working good
        # Same idea but with boolean vars
        done = True
        for i in range(len(self.vehicles)):
            v = self.vehicles[i].p.v
            e = self.epsilon
            # if they are all equal or less than
            condition = (v <= e or isclose(v, e, abs_tol=1e-5))
            # if(condition == False):
            #     print(v,e)

            done = done and condition
        return done

    def get_mm1(self, i):
        l = len(self.vehicles)
        m, m1 = self.vehicles[i%l], self.vehicles[(i+1)%l]
        pm, pm1 = m.p, m1.p
        return pm, pm1

    def _compute_dt(self):
        dt = 99e99
        a_max = -99e99
        # Get the maximum of the acceleration
        for i in range(len(self.vehicles)):
            a_max = max(a_max, self.vehicles[i].p.a)
        for i in range(len(self.vehicles)):
            pm, pm1 = self.get_mm1(i)
            vm = self.vehicles[i].p.v#pm.v
            dxm = Vehicle.compute_dx(pm, pm1, len(self.vehicles), self.track_length)
            if(vm > 0):
                d = dxm / vm
                if(a_max > 0):
                    val = (-vm + sqrt(vm*vm+2*dxm*a_max))/(a_max)
                    dt = min(dt, val)
                else:
                    dt = min(d,dt)
        if((self.timesteps // self.tau) < (self.timesteps+dt) // self.tau):            
            dt = self.tau * (self.timesteps//self.tau+1) - self.timesteps
        
        return dt        

    def _compute_entrypoint(self):
        x0 = self.vehicles[0].p.x   # zero index is the latest car
        xM = self.vehicles[-1].p.x
        entrypoint = 0.0
        # If we have more than one car and the leader car crossed the origin (start line)
        #   Then the leader car will have less value than the latest one
        #   In order to calculate the mean correctly
        if(x0 >= xM and len(self.vehicles) > 1):
            entrypoint = xM + (x0 - xM)/2 
        else:
            entrypoint = xM + (x0 - xM + self.track_length) / 2

        if(entrypoint > self.track_length):
            entrypoint -= self.track_length
        return entrypoint

    def _add_vehicle(self):
        p0 = self.vehicles[0].p   # zero index is the latest car
        pM = self.vehicles[-1].p
        dx = Vehicle.compute_dx(p0, pM, len(self.vehicles), self.track_length)
        # After more than two vehicles are moving on the road a new vehicle appears at the point with the coordinate midway between 
        # 𝑥1(𝑡) 𝑎𝑛𝑑 𝑥𝑀(𝑡) vehicles if the distance between them: |𝑥1(𝑡) − 𝑥𝑀(𝑡)| ≥ 7 [𝑚]. 
        if(self.timesteps // self.tau < (self.timesteps+self.dt) // self.tau and abs(dx) >= 7):
            self.id -= 1
            entrypoint = self._compute_entrypoint()
            v = Vehicle(x0=entrypoint, v0=self.initial_velocity, a0=self.initial_acceleration, id=self.id)
            # Add a vehicle to the begining of the list in order to have the correc indexing
            # indexing [0,1,2,3] where 3 is the oldest vehicles
            self.vehicles.insert(0, v)
            self.add_flag = True

    def _update_vehicles(self):
        for i in range(len(self.vehicles)):
            dt = self.dt
            self.vehicles[i].update(dt, self.epsilon, self.track_length)
    
    def _update_vehicles_acc(self):
        for i in range(len(self.vehicles)):
            pm, pm1 = self.get_mm1(i)
            self.vehicles[i].update_acc(pm, pm1, len(self.vehicles), self.track_length)

    def _print_report(self, sleep_time=0.01):
        print(f"dt: {self.dt}\ttimesteps: {self.timesteps}\tNum. Steps: {self.num_steps}")
        print(f"Num. Vehicles: {len(self.vehicles)}")
        print("------------------")
        sleep(sleep_time)
        
    def _speed_limit(self):
        # We speed limit all the followers if we added on in the middle
        flag = -1*(self.add_flag==True)
        self.add_flag = False
        for i in range(len(self.vehicles)+flag):
            pm, pm1 = self.get_mm1(i)
            dx = Vehicle.compute_dx(pm, pm1, len(self.vehicles), track_length)
            # If ∀𝑚 Δ𝑥𝑚(𝑡) < 7 [м], then 𝑣𝑚(𝑡) = 0
            if(dx < 7 or self.vehicles[i].p.v < 0.0):
                self.vehicles[i].p.v = 0.0

    def _update_counters(self):
        self.timesteps += self.dt
        self.num_steps += 1

    def update(self):
        # To compute dt we need to update acceleration first
        self._update_vehicles_acc()
        self.dt = self._compute_dt()
        self._update_vehicles()
        self._add_vehicle()
        self._speed_limit()
        done = self._compute_termination()
        self._update_counters()
        # self._print_report()
        return done, self.timesteps, self.vehicles

def print_report(data, file="output.txt"):
    timesteps = data
    with open(file, 'w+') as f:
        f.write(str(timesteps))

def parse_file(file="input.txt", line_num=0):
    with open(file, 'r') as f:
        lines = f.readlines()
        data = list(lines[line_num].split(' ')) # only one line
        data = [float(d) for d in data]
    return data # track_length, tau, epsilon

if __name__ == '__main__':
    track_length, tau, epsilon = parse_file()
    traffic = Traffic(
                        track_length=track_length,
                        tau=tau,
                     )

    traffic.reset()
    done = False
    steps = 0
    while not done:
        done, timesteps, vehicles = traffic.update()
        steps += 1
    print_report(timesteps)
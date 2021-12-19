# Microscopic model of traffic flow
# Sources:
# - https://gafferongames.com/post/integration_basics/
# - https://towardsdatascience.com/simulating-traffic-flow-in-python-ee1eab4dd20f
# Run it: python3 main.py

# TODO: take into consideration the possibility of having only 1 vehicle
from math import sqrt
from time import sleep

class Physics:
    def __init__(self, x0=0.0, v0=0.0, a0=0.0):
        self.x = x0    # starting point
        self.v = v0
        self.a = a0


class Vehicle:
    def __init__(self, x0=0.0, v0=0.0, a0=0.0, id=None):
        self.p = Physics(x0, v0, a0)
        self.id = id

    def _compute_c_coeff(self, pm, pm1):
        xm, xm1 = pm.x, pm1.x
        dxm = xm1 - xm
        c = 0.0
        if(dxm >= 45):  # Free-flow
            c = 1075/dxm
        elif(dxm >= 25):   # Synchronized-flow
            c = 0.6*dxm + 362/dxm
        elif(dxm >= 7):
            c = 0.6*dxm
        return c

    def _compute_dx(self, pm, pm1, len_traffic, track_length):
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
    
    def update(self, pm, pm1, dt, len_traffic, eps, track_length):
        dv = pm1.v - pm.v
        dx = self._compute_dx(pm, pm1, len_traffic, track_length)
        old_x = self.p.x
        c = self._compute_c_coeff(pm, pm1)
        self.p.a = c*dv/dx

        # Semi-implicit Euler integration: https://www.wikiwand.com/en/Semi-implicit_Euler_method
        # "Since this is only an approximation, the speed can become negative at times (but the model does not allow for that). An instability arises when the speed is negative, and the position and speed diverge into negative infinity"
        # Source: https://towardsdatascience.com/simulating-traffic-flow-in-python-ee1eab4dd20f
        # "To overcome this problem, whenever we predict a negative speed we will set it equal to zero and work out way from there:"
        if(self.p.v + self.p.a*dt < 0):
            self.p.x -= 0.5*self.p.v*self.p.v/self.p.a
            self.p.v = 0.0
        else:
            self.p.v += self.p.a*dt
            self.p.x += self.p.v*dt + self.p.a*dt*dt/2


        # If the new position is nearly equal
        if(self.p.x - (dx+old_x) < eps and dx+old_x < self.p.x):
            self.p.x = (dx+old_x)
        
        # Ring round track -> limit the position
        if(self.p.x >= track_length):
            self.p.x = self.p.x - self.p.x // track_length * track_length

    def get_attributes(self):
        return self.p, self.id

# Discrete microscopic traffic modeling
class Traffic:
    def __init__(self, track_length=None, tau=5, initial_velocity=40, initial_acceleration=0.0, epsilon=0.0001):
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
        self.vehicles = [Vehicle(v0=self.initial_velocity, a0=self.initial_acceleration, id=self.id)]
        self.timesteps = self.tau # As we added a vehcile, and as initially t=0 the road is empty
        self.entrypoint = 0.0
        self.last_time_add = 0.0
        self.num_steps = 0

    def _compute_termination(self):
        done = True
        for i in range(len(self.vehicles)):
            condition = self.vehicles[i].p.v <= self.epsilon
            done = done and condition
        return done

    # Get attributes of m and m1 vechiles
    def get_mm1(self, m):
        l = len(self.vehicles)
        m, m1 = self.vehicles[m%l], self.vehicles[(m+1)%l]
        # m, m1 = self.vehicles[m], self.vehicles[m-1]    # As the zeroth index is the leading car and the car with bigger index is the follower
        pm, pm1 = m.get_attributes()[0], m1.get_attributes()[0]
        return pm, pm1

    def _compute_dt(self):
        dt = 99e99
        a_max = -99e99
        for i in range(len(self.vehicles)):
            a_max = max(a_max, self.vehicles[i].p.a)
        for i in range(len(self.vehicles)):
            pm, pm1 = self.get_mm1(i)
            vm = pm.v
            dxm = pm1.x - pm.x
            if(vm > 0):
                d = dxm / vm
                if(a_max > 0):
                    val = (-vm + sqrt(vm*vm+4*dxm*a_max))/(2*a_max)
                    dt = min(dt, val)
                else:
                    dt = min(d,dt)
        return dt        
        # if(len(self.vehicles) > 1):
        #     dt = 99e99
        #     a_max = -99e99
        #     for i in range(len(self.vehicles)):
        #         a_max = max(a_max, self.vehicles[i].p.a)
        #     for i in range(len(self.vehicles)):
        #         pm, pm1 = self.get_mm1(i)
        #         vm = pm.v
        #         dxm = pm1.x - pm.x
        #         val = (-vm + sqrt(vm*vm+4*dxm*a_max))/(2*a_max)
        #         dt = min(dt, val)
        #     self.dt = dt

    def _compute_entrypoint(self):
        p1 = self.vehicles[-1]
        pM = self.vehicles[0]
        # TODO: Check this condition 
        if(len(self.vehicles) > 2 and abs(p1.x - pM.x) >= 7):
            p1 = self.vehicles[-1]
            pM = self.vehicles[0]
            self.entrypoint = (pM.x + p1.x)/2

    def _add_vehicle(self):
        self._compute_entrypoint()
        if(self.last_time_add > self.tau):
            self.id -= 1
            v = Vehicle(x0=self.entrypoint, v0=self.initial_velocity, a0=self.initial_acceleration, id=self.id)
            # self.vehicles.append(v)
            # Add vehicle to the begining of the list in order to have the correc indexing
            # indexing [0,1,2,3] where 3 is the oldest vehicles
            self.vehicles.insert(0, v)
            self.last_time_add = 0.0

    def _update_vehicles(self):
        for i in range(len(self.vehicles)):
            pm, pm1 = self.get_mm1(i)
            dt = self.dt
            self.vehicles[i].update(pm, pm1, dt, len(self.vehicles), self.epsilon, self.track_length)

    def _print_report(self, sleep_time=0.1):
        print(f"dt: {self.dt}\ttimesteps: {self.timesteps}\tNum. Steps: {self.num_steps}")
        print(f"Num. Vehicles: {len(self.vehicles)}")
        sleep(sleep_time)

    def update(self):
        self.dt = self._compute_dt()
        self._update_vehicles()
        self._add_vehicle()
        done = self._compute_termination()
        self._print_report()
        self.timesteps += self.dt
        self.num_steps += 1
        self.last_time_add += self.dt

        return done, self.timesteps, self.vehicles

def print_report(data, file="output.txt"):
    timesteps, correct_timesteps, vehicles = data
    error = abs(timesteps - correct_timesteps) / max(1, abs(timesteps))
    print(f"Timesteps: {timesteps},\tCorrect: {correct_timesteps},\tError: {error}\nNum. vehicles: {len(vehicles)}")
    with open(file, 'w+') as f:
        f.write(str(timesteps))

def parse_file(file="input.txt", line_num=0):
    with open(file, 'r') as f:
        lines = f.readlines()
        data = list(lines[line_num].split(' ')) # only one line
        # data = map(lambda x: float(x), data)
        data = [float(d) for d in data]
    return data # track_length, tau, epsilon

if __name__ == '__main__':
    track_length, tau, epsilon = parse_file()
    # print(track_length, tau, epsilon)
    traffic = Traffic(
                        track_length=track_length,
                        tau=tau,
                        # initial_velocity=initial_velocity,
                        # initial_acceleration=initial_acceleration
                     )

    traffic.reset()
    done = False
    while not done:
        done, timesteps, vehicles = traffic.update()

    # correct_timesteps = timesteps
    correct_timesteps, num_vehicles = parse_file(file="correct.txt")
    # print(correct_timesteps)
    print(f"Correct number of vehicles: {num_vehicles}")
    print_report([timesteps, correct_timesteps, vehicles])
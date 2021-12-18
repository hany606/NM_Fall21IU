class Vehicle:
    def __init__(self, x0=0.0, v0=0.0, a0=0.0, id=0):
        self.x = x0    # starting point
        self.v = v0
        self.a = a0
        self.id = id

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
        self.vehicles = []
        self.timesteps = 0.0
        self.dt = 0.0
        
    def _compute_c_coeff(self, xm1, xm):
        dxm = xm1 - xm
        c = None
        if(dxm >= 45):  # Free-flow
            c = 1075/dxm
        elif(25 <= dxm < 45):   # Synchronized-flow
            c = 0.6*dxm + 362/dxm
        elif(7 <= dxm < 25):
            c = 0.6*dxm
        

    def _compute_termination(self):
        return True

    def update(self):
        self.timesteps += self.dt
        done = self._compute_termination()

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
    correct_timesteps = parse_file(file="correct.txt")[0]
    # print(correct_timesteps)
    print_report([timesteps, correct_timesteps, vehicles])
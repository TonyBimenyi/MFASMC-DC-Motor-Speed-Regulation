import numpy as np
import random
import sys
import matplotlib.pyplot as plt

# Simple MFAC + SMC for Chatbot Response Adaptation
# - MFAC adapts 'relevance' estimate (phi) based on user feedback.
# - SMC ensures robust 'satisfaction' tracking to target (0.8).
# - Chatbot: Rule-based responses; control adjusts 'confidence' to select better ones.

class MFAC_SMC_Chatbot:
    def __init__(self):
        self.phi = 1.0  # Initial pseudo-gradient (relevance sensitivity)
        self.eta = 0.5  # MFAC adaptation step
        self.rho = 1.0  # MFAC control gain
        self.mu = 1e-3  # Regularization
        self.lambda_ = 1e-4
        self.alpha = 1.0  # SMC surface slope
        self.epsilon = 0.5  # SMC reaching gain
        self.beta = 0.1  # SMC weighting
        self.target_satisfaction = 0.8  # Target user satisfaction
        self.satisfaction = 0.5  # Initial
        self.e_prev = 0  # Previous error for SMC
        self.satisfaction_prev = 0.5
        self.confidence = 0.5  # Initial confidence for response selection
        self.confidence_prev = 0.5
        self.responses = {
            'hello': ['Hi there!', 'Hello! How can I help?', 'Hey! What\'s up?'],
            'how are you': ['I\'m good, thanks!', 'Doing great! And you?'],
            'bye': ['Goodbye!', 'See you later!', 'Bye! Have a good day.'],
            'default': ['Interesting, tell me more.', 'I see. What else?', 'That makes sense.']
        }
        self.satisfaction_history = [self.satisfaction]  # For plotting
        self.confidence_history = [self.confidence]     # For plotting
        self.phi_history = [self.phi]                   # For plotting
        self.turns = [0]                                # Turn counter for x-axis
    
    def get_response(self, user_input):
        # Simple rule-based response (confidence could bias selection in advanced version)
        user_input_lower = user_input.lower()
        for key in self.responses:
            if key in user_input_lower:
                return random.choice(self.responses[key])
        return random.choice(self.responses['default'])
    
    def adapt_with_mfac_smc(self, user_feedback):  # Feedback: 1 (good), 0 (bad)
        # Update satisfaction (simple moving average)
        self.satisfaction = 0.7 * self.satisfaction + 0.3 * user_feedback
        e = self.target_satisfaction - self.satisfaction  # Error
        # MFAC part: Adapt phi and compute delta_conf_mfac
        du = self.confidence - self.confidence_prev
        d_satisfaction = self.satisfaction - self.satisfaction_prev
        self.phi += self.eta * (d_satisfaction - self.phi * du) / (self.mu + du**2)
        denom = self.lambda_ + self.phi**2
        delta_conf_mfac = (self.rho * self.phi * e) / denom
        # SMC part: Sliding surface and robust correction
        e_dot = (e - self.e_prev) / 1.0  # Approx dt=1 (per turn)
        s = e + self.alpha * self.e_prev  # Surface
        y_dot = -0.01 * self.satisfaction  # Approx satisfaction derivative
        yd_dot = 0  # Target constant
        equiv_u = (yd_dot - y_dot + self.alpha * e_dot) / self.phi
        reach_u = self.epsilon * np.sign(s)
        u_smc = (equiv_u + reach_u) / self.phi
        # Total adjustment to confidence
        delta_conf_total = delta_conf_mfac + self.beta * u_smc
        self.confidence += delta_conf_total * 0.1  # Scaled for stability
        self.confidence = np.clip(self.confidence, 0, 1)
        # Store previous values
        self.satisfaction_prev = self.satisfaction
        self.e_prev = e
        self.confidence_prev = self.confidence
        
        # Log for plotting (after adaptation)
        self.turns.append(len(self.turns))
        self.satisfaction_history.append(self.satisfaction)
        self.confidence_history.append(self.confidence)
        self.phi_history.append(self.phi)
    
    def chat(self, user_input, feedback=None):
        response = self.get_response(user_input)
        if feedback is not None:
            self.adapt_with_mfac_smc(feedback)
        return response
    
    def plot_adaptation(self):
        # Plot evolution over turns
        fig, axs = plt.subplots(3, 1, figsize=(10, 8))
        
        axs[0].plot(self.turns, self.satisfaction_history, 'g-', linewidth=2, label='Satisfaction')
        axs[0].axhline(y=self.target_satisfaction, color='r', linestyle='--', label='Target')
        axs[0].set_ylabel('Satisfaction')
        axs[0].legend()
        axs[0].grid(True)
        
        axs[1].plot(self.turns, self.confidence_history, 'b-', linewidth=2, label='Confidence')
        axs[1].set_ylabel('Confidence')
        axs[1].legend()
        axs[1].grid(True)
        
        axs[2].plot(self.turns, self.phi_history, 'm-', linewidth=2, label='Phi (Adapted)')
        axs[2].set_ylabel('Phi')
        axs[2].set_xlabel('Turns')
        axs[2].legend()
        axs[2].grid(True)
        
        plt.suptitle('MFAC + SMC Adaptation Over Chat Turns')
        plt.tight_layout()
        plt.show()

# Interactive Demo
if __name__ == "__main__":
    bot = MFAC_SMC_Chatbot()
    # print("=== MFAC + SMC Adaptive Chatbot Demo ===")
    # print("Type messages. After each response, rate 1 (good) or 0 (bad). Type 'quit' to stop.")
    # print("Bot adapts satisfaction/confidence over time. Plots at end.\n")
    
    # while True:
    #     try:
    #         user = input("You: ").strip()
    #         if user.lower() == 'quit':
    #             break
    #         resp = bot.chat(user)
    #         print(f"Bot: {resp}")
    #         fb_input = input("Feedback (1 good/0 bad): ").strip()
    #         if fb_input in ['1', '0']:
    #             fb = int(fb_input)
    #             bot.chat(user, fb)  # Trigger adaptation (dummy call)
    #             print(f"Updated: Satisfaction={bot.satisfaction:.2f}, Confidence={bot.confidence:.2f}, Phi={bot.phi:.2f}\n")
    #         else:
    #             print("Invalid feedback. Skipping adaptation.\n")
    #     except KeyboardInterrupt:
    #         break
    
    # print(f"\nFinal Stats: Satisfaction={bot.satisfaction:.2f}, Confidence={bot.confidence:.2f}")
    # print("Phi (adapted sensitivity):", bot.phi)
    bot.plot_adaptation()  # Generate and show plots
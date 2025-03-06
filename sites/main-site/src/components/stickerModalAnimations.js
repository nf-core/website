/**
 * Sticker Modal Animations
 * Provides enhanced animations for the sticker modal
 */

// Animation for opening the modal with a sticker
export function animateModalOpen(modal, modalContainer, modalHex) {
  if (!modal || !modalContainer || !modalHex) return;

  // Reset any existing animations
  modalContainer.style.animation = '';
  modalHex.style.animation = '';

  // Apply opening animations
  modalContainer.style.animation = 'modalContainerIn 0.5s cubic-bezier(0.175, 0.885, 0.32, 1.275) forwards';
  modalHex.style.animation = 'modalHexIn 0.7s cubic-bezier(0.175, 0.885, 0.32, 1.275) forwards';
}

// Animation for closing the modal
export function animateModalClose(modal, modalContainer, modalHex, callback) {
  if (!modal || !modalContainer || !modalHex) {
    if (callback) callback();
    return;
  }

  // Apply closing animations
  modalContainer.style.animation = 'modalContainerOut 0.4s cubic-bezier(0.55, 0.085, 0.68, 0.53) forwards';
  modalHex.style.animation = 'modalHexOut 0.3s cubic-bezier(0.55, 0.085, 0.68, 0.53) forwards';

  // Execute callback after animation completes
  setTimeout(() => {
    if (callback) callback();
  }, 400);
}

// Add confetti effect when opening a sticker
export function addConfettiEffect(container) {
  if (!container) return;

  // Create confetti elements
  const confettiCount = 50;
  const colors = ['#f94144', '#f3722c', '#f8961e', '#f9c74f', '#90be6d', '#43aa8b', '#577590'];

  // Remove any existing confetti
  const existingConfetti = container.querySelectorAll('.confetti');
  existingConfetti.forEach(c => c.remove());

  // Create new confetti pieces
  for (let i = 0; i < confettiCount; i++) {
    const confetti = document.createElement('div');
    confetti.className = 'confetti';
    confetti.style.backgroundColor = colors[Math.floor(Math.random() * colors.length)];
    confetti.style.left = Math.random() * 100 + '%';
    confetti.style.top = (Math.random() * 20 - 10) + '%';
    confetti.style.transform = `rotate(${Math.random() * 360}deg)`;
    confetti.style.width = (Math.random() * 8 + 4) + 'px';
    confetti.style.height = (Math.random() * 6 + 2) + 'px';
    confetti.style.opacity = Math.random() + 0.5;

    // Random animation duration and delay
    const duration = Math.random() * 3 + 2;
    const delay = Math.random() * 0.5;
    confetti.style.animation = `confettiFall ${duration}s ${delay}s ease-in forwards`;

    container.appendChild(confetti);

    // Remove confetti after animation completes
    setTimeout(() => {
      confetti.remove();
    }, (duration + delay) * 1000);
  }
}

// Add a shine effect to the sticker
export function addShineEffect(hexElement) {
  if (!hexElement) return;

  // Create shine element if it doesn't exist
  let shine = hexElement.querySelector('.shine-effect');
  if (!shine) {
    shine = document.createElement('div');
    shine.className = 'shine-effect';
    hexElement.appendChild(shine);
  }

  // Reset animation
  shine.style.animation = 'none';

  // Trigger reflow
  void shine.offsetWidth;

  // Start animation
  shine.style.animation = 'shineEffect 2s ease-in-out';
}

// Add a pulsing effect to highlight the sticker
export function addPulseEffect(hexElement) {
  if (!hexElement) return;

  // Add pulse class
  hexElement.classList.add('pulse-effect');

  // Remove class after animation completes
  setTimeout(() => {
    hexElement.classList.remove('pulse-effect');
  }, 1000);
}

// Add CSS styles for the animations
export function injectAnimationStyles() {
  // Check if styles are already added
  if (document.getElementById('sticker-animation-styles')) return;

  const styleElement = document.createElement('style');
  styleElement.id = 'sticker-animation-styles';

  styleElement.textContent = `
    @keyframes modalContainerIn {
      0% { transform: scale(0.8); opacity: 0; }
      100% { transform: scale(1); opacity: 1; }
    }

    @keyframes modalContainerOut {
      0% { transform: scale(1); opacity: 1; }
      100% { transform: scale(0.8); opacity: 0; }
    }

    @keyframes modalHexIn {
      0% { transform: scale(0.5) rotate(-15deg); opacity: 0; }
      40% { transform: scale(1.1) rotate(5deg); opacity: 1; }
      70% { transform: scale(0.95) rotate(-2deg); }
      100% { transform: scale(1) rotate(0); }
    }

    @keyframes modalHexOut {
      0% { transform: scale(1) rotate(0); opacity: 1; }
      100% { transform: scale(0.5) rotate(15deg); opacity: 0; }
    }

    @keyframes confettiFall {
      0% { transform: translateY(-10px) rotate(0deg); opacity: 1; }
      100% { transform: translateY(500px) rotate(360deg); opacity: 0; }
    }

    @keyframes shineEffect {
      0% { left: -100%; opacity: 0; }
      20% { opacity: 0.8; }
      100% { left: 100%; opacity: 0; }
    }

    .confetti {
      position: absolute;
      z-index: 100;
      pointer-events: none;
    }

    .shine-effect {
      position: absolute;
      top: 0;
      left: -100%;
      width: 50%;
      height: 100%;
      background: linear-gradient(
        to right,
        rgba(255, 255, 255, 0) 0%,
        rgba(255, 255, 255, 0.8) 50%,
        rgba(255, 255, 255, 0) 100%
      );
      transform: skewX(-20deg);
      pointer-events: none;
      z-index: 10;
    }

    .pulse-effect {
      animation: pulse 1s cubic-bezier(0.175, 0.885, 0.32, 1.275);
    }

    @keyframes pulse {
      0% { transform: scale(1); }
      50% { transform: scale(1.1); }
      100% { transform: scale(1); }
    }
  `;

  document.head.appendChild(styleElement);
}

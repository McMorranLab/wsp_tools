var currentSlide = 0;
var slides = document.getElementsByClassName("slide");
var currentAnimation = [];
for (var i = 0; i < slides.length; i++) {
	currentAnimation.push(0);
	slides[i].classList.add("hide");
	var slideNum = document.createElement("P");
	slideNum.innerHTML = String(i);
	slideNum.style = "position: absolute; bottom: 20px; right: 40px; font-size:20px;"
	slides[i].appendChild(slideNum);
}
slides[0].classList.remove("hide");

function lastSlide() {
	if (currentSlide + 1 < slides.length) {
		return (false);
	}
	return (true);
}

function firstSlide() {
	if (currentSlide > 0) {
		return (false);
	}
	return (true);
}

function hasAnimations() {
	if (slides[currentSlide].querySelectorAll(".animate-go, .animate-stay").length > 0) {
		return (true);
	} else {
		return (false);
	}
}

function lastAnimation() {
	var animations = slides[currentSlide].querySelectorAll(".animate-go, .animate-stay");
	for (var i = 0; i < animations.length; i++) {
		if (animations[i].classList.contains(String(currentAnimation[currentSlide] + 1))) {
			return (false);
		}
	}
	return (true);
}

function firstAnimation() {
	if (currentAnimation[currentSlide] > 0) {
		return (false);
	}
	return (true);
}

function increaseSlide() {
	slides[currentSlide].classList.add("hide");
	currentSlide += 1;
	slides[currentSlide].classList.remove("hide");
	var animations = slides[currentSlide].querySelectorAll(".animate-go, .animate-stay");
	for (var i = 0; i < animations.length; i++) {
		if (!animations[i].classList.contains("0")) {
			animations[i].classList.add("hide");
		}
	}
}

function decreaseSlide() {
	slides[currentSlide].classList.add("hide");
	currentSlide -= 1;
	slides[currentSlide].classList.remove("hide");
}

function decreaseAnimation() {
	var animateGo = slides[currentSlide].getElementsByClassName("animate-go " + String(currentAnimation[currentSlide]));
	var animateStay = slides[currentSlide].getElementsByClassName("animate-stay " + String(currentAnimation[currentSlide]));
	for (var i = 0; i < animateGo.length; i++) {
		animateGo[i].classList.add("hide");
	}
	for (var i = 0; i < animateStay.length; i++) {
		animateStay[i].classList.add("hide");
	}
	currentAnimation[currentSlide] -= 1;
	animateGo = slides[currentSlide].getElementsByClassName("animate-go " + String(currentAnimation[currentSlide]));
	for (var i = 0; i < animateGo.length; i++) {
		animateGo[i].classList.remove("hide");
	}
}

function increaseAnimation() {
	var animateGo = slides[currentSlide].getElementsByClassName("animate-go " + String(currentAnimation[currentSlide]));
	for (var i = 0; i < animateGo.length; i++) {
		animateGo[i].classList.add("hide");
	}
	currentAnimation[currentSlide] += 1;
	animateGo = slides[currentSlide].getElementsByClassName("animate-go " + String(currentAnimation[currentSlide]));
	var animateStay = slides[currentSlide].getElementsByClassName("animate-stay " + String(currentAnimation[currentSlide]));
	for (var i = 0; i < animateGo.length; i++) {
		animateGo[i].classList.remove("hide");
	}
	for (var i = 0; i < animateStay.length; i++) {
		animateStay[i].classList.remove("hide");
	}
}

window.addEventListener("keydown", function(e) {
	if ((e.keyCode === 39) || (e.keyCode === 32)) {
		if (hasAnimations()) {
			if (!lastAnimation()) {
				increaseAnimation();
			} else if (!lastSlide()) {
				increaseSlide();
			}
		} else if (!lastSlide()) {
			increaseSlide();
		}
	} else if (e.keyCode === 37) {
		if (hasAnimations()) {
			if (!firstAnimation()) {
				decreaseAnimation();
			} else if (!firstSlide()) {
				decreaseSlide();
			}
		} else if (!firstSlide()) {
			decreaseSlide();
		}
	}
})
